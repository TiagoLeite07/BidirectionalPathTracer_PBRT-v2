#include "stdafx.h"
#include "integrators/bidirectionalpath.h"
#include "scene.h"
#include "intersection.h"
#include "paramset.h"
#include "camera.h"
#include "montecarlo.h"
#include "progressreporter.h"

void BidirectionalPathIntegrator::RequestSamples(Sampler *sampler, Sample *sample, 
												 const Scene *scene) {
	for (int i = 0; i < SAMPLE_DEPTH; ++i) {
        lightSampleOffsets[i] = LightSampleOffsets(1, sample);
        bsdfSampleOffsets[i] = BSDFSampleOffsets(1, sample);
        pathSampleOffsets[i] = BSDFSampleOffsets(1, sample);
		lightpathSampleOffsets[i] = BSDFSampleOffsets(1, sample);
    }
}

Spectrum BidirectionalPathIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &r, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
	// Declare common path integration variables
    Spectrum pathThroughput = 1., L = 0.;
    RayDifferential ray(r);
    bool specularBounce = false;
    Intersection localIsect;
    const Intersection *isectp = &isect;
	int lp;
	vector<LightPathPoint> lightpath;

	//-------LIGHT PATH-------

	if (scene->lights.size() == 0) return 0;
	// Compute samples for emitted rays from lights
	vector<float> lightNum(1);
	vector<float> lightSampPos(2, 0.f);
	vector<float> lightSampComp(1, 0.f);
	vector<float> lightSampDir(2, 0.f);
	LDShuffleScrambled1D(1, 1, &lightNum[0], rng);
	LDShuffleScrambled2D(1, 1, &lightSampPos[0], rng);
	LDShuffleScrambled1D(1, 1, &lightSampComp[0], rng);
	LDShuffleScrambled2D(1, 1, &lightSampDir[0], rng);

	// Precompute information for light sampling densities
	Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);

	int sampOffset = 0;

	// Choose light source to trace virtual light path from
	float lightPdf;
	int ln = lightDistribution->SampleDiscrete(lightNum[sampOffset],
													&lightPdf);
	Light *light = scene->lights[ln];

	// Sample ray leaving light source for light path
	RayDifferential lightray;
	float pdf;
	LightSample ls(lightSampPos[2*sampOffset], lightSampPos[2*sampOffset+1],
					   lightSampComp[sampOffset]);
	Normal Nl;
	Spectrum alpha = light->Sample_L(scene, ls, lightSampDir[2*sampOffset],
									 lightSampDir[2*sampOffset+1],
									 0, &lightray, &Nl, &pdf); //camera->shutterOpen=0
	if (!(pdf == 0.f || alpha.IsBlack())) {
		alpha /= pdf * lightPdf;
		Intersection lightisect;
		int lightbounces = 0;
		while (scene->Intersect(lightray, &lightisect) && !alpha.IsBlack()) {
			// Create light path point and sample new ray for path
			alpha *= renderer->Transmittance(scene, RayDifferential(lightray), NULL,
											 rng, arena);
			Vector wo = -lightray.d;
			BSDF *bsdf = lightisect.GetBSDF(lightray, arena);

			// Create light path point at ray intersection point
			Spectrum contrib = alpha * bsdf->rho(wo, rng) / M_PI;
			lightpath.push_back(LightPathPoint(lightisect.dg.p, lightisect.dg.nn, contrib,
													lightisect.rayEpsilon));

			// Sample new ray direction and update weight for light path
			Vector wi;
			float pdf;
			BSDFSample bsdfSample;
			if (lightbounces < SAMPLE_DEPTH)
				bsdfSample = BSDFSample(sample, lightpathSampleOffsets[lightbounces], 0);
			else
				bsdfSample = BSDFSample(rng);
			Spectrum fr = bsdf->Sample_f(wo, &wi, bsdfSample, &pdf);
			if (fr.IsBlack() || pdf == 0.f)
				break;
			Spectrum contribScale = fr * AbsDot(wi, bsdf->dgShading.nn) / pdf;

			// Possibly terminate light path with Russian roulette
			float rrProb = min(1.f, contribScale.y());
			if (rng.RandomFloat() > rrProb)
				break;
			if (lightbounces == maxDepth)
				break;
			alpha *= contribScale / rrProb;
			lightray = RayDifferential(lightisect.dg.p, wi, lightray, lightisect.rayEpsilon);
			lightbounces++;
		}
		arena.FreeAll();
	}
	delete lightDistribution;

	//------CAMERA PATH-------

	for (int bounces = 0; bounces < SAMPLE_DEPTH ; ++bounces) {
        // Possibly add emitted light at path vertex
        if (bounces == 0 || specularBounce)
            L += pathThroughput * isectp->Le(-ray.d);

        // Sample illumination from lights to find path contribution
        BSDF *bsdf = isectp->GetBSDF(ray, arena);
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;
        Vector wo = -ray.d;

		if(bounces < SAMPLE_DEPTH){
			L += pathThroughput * UniformSampleAllLights(scene, renderer, arena, p, n,
						wo, isectp->rayEpsilon, ray.time, bsdf, sample, rng,
						lightSampleOffsets, bsdfSampleOffsets);
		}

		Spectrum L_bdpt = 0.;
		//if (!specularBounce) {
        for (uint32_t i = 0; i < lightpath.size(); ++i) {
			const LightPathPoint &lpp = lightpath[i];
			// Compute light's path tentative contribution _Llight_
			float d2 = DistanceSquared(p, lpp.p);
			Vector wi_bdpt = Normalize(lpp.p - p);
			float G = AbsDot(wi_bdpt, n) * AbsDot(wi_bdpt, lpp.n) / d2;
			G = min(G, gLimit);
			Spectrum f = bsdf->f(wo, wi_bdpt);
			if (G == 0.f || f.IsBlack()) continue;
			Spectrum Llight = f * G * lpp.pathContrib;
			RayDifferential connectRay(p, wi_bdpt, ray, isect.rayEpsilon,
										   sqrtf(d2) * (1.f - lpp.rayEpsilon));
			Llight *= renderer->Transmittance(scene, connectRay, NULL, rng, arena);

			// Add contribution from _LightPathPoint_ _lpp_ to _L_bdpt_
			if (!scene->IntersectP(connectRay)) {
				L_bdpt += pathThroughput * Llight;
			}
		}
		//L += L_bdpt;
		// Add mean contribution of _L_bdpt_ to _L_
		lp = lightpath.size();
		if (lp != 0)
			L += (L_bdpt/lp);
		//}
		//else {
			//lp=0;
		//}
		//lp=0;



        // Sample BSDF to get new path direction

        // Get _outgoingBSDFSample_ for sampling new path direction
        BSDFSample outgoingBSDFSample;
        if (bounces < SAMPLE_DEPTH)
            outgoingBSDFSample = BSDFSample(sample, pathSampleOffsets[bounces],
                                            0);
        else
            outgoingBSDFSample = BSDFSample(rng);
        Vector wi;
        float pdf;
        BxDFType flags;
        Spectrum f = bsdf->Sample_f(wo, &wi, outgoingBSDFSample, &pdf,
                                    BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == 0.)
            break;
        specularBounce = (flags & BSDF_SPECULAR) != 0;
		if(specularBounce)
			pathThroughput *= f * AbsDot(wi, n) / (pdf);
		else
			pathThroughput *= f * AbsDot(wi, n) / (pdf * (lp+1));
        ray = RayDifferential(p, wi, ray, isectp->rayEpsilon);

        // Possibly terminate the path
        if (bounces > 3) {
            float continueProbability = min(.5f, pathThroughput.y());
            if (rng.RandomFloat() > continueProbability)
                break;
            pathThroughput /= continueProbability;
        }
        if (bounces == maxDepth)
            break;

        // Find next vertex of path
        if (!scene->Intersect(ray, &localIsect)) {
            if (specularBounce)
                for (uint32_t i = 0; i < scene->lights.size(); ++i)
                   L += pathThroughput * scene->lights[i]->Le(ray);
            break;
        }
        pathThroughput *= renderer->Transmittance(scene, ray, NULL, rng, arena);
        isectp = &localIsect;
    }
	lightpath.clear();
    return L;
}

BidirectionalPathIntegrator *CreateBidirectionalPathSurfaceIntegrator(const ParamSet &params) {
	int maxDepth = params.FindOneInt("maxdepth", 5);
	float glimit = params.FindOneFloat("glimit", 10.f);
    return new BidirectionalPathIntegrator(maxDepth, glimit);
}