#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_BIDIRECTIONALPATH_H
#define PBRT_INTEGRATORS_BIDIRECTIONALPATH_H

#include "pbrt.h"
#include "integrator.h"

struct LightPathPoint {
    LightPathPoint() { }
    LightPathPoint(const Point &pp, const Normal &nn, const Spectrum &c,
                 float reps)
        : p(pp), n(nn), pathContrib(c), rayEpsilon(reps) { }
    Point p;
    Normal n;
    Spectrum pathContrib;
    float rayEpsilon;
};

class BidirectionalPathIntegrator : public SurfaceIntegrator {
public:
	Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;

	void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);

	//vector<LightPathPoint> GenerateLightPath(const Scene *scene, const Renderer *renderer);

	BidirectionalPathIntegrator(int md, float gl) { 
		maxDepth = md;
		gLimit = gl;
	}

private:
	int maxDepth;
	float gLimit;
#define SAMPLE_DEPTH 5
    LightSampleOffsets lightSampleOffsets[SAMPLE_DEPTH];
    BSDFSampleOffsets bsdfSampleOffsets[SAMPLE_DEPTH];
    BSDFSampleOffsets pathSampleOffsets[SAMPLE_DEPTH];
	BSDFSampleOffsets lightpathSampleOffsets[SAMPLE_DEPTH];

	//vector<LightPathPoint> lightpath;
};

BidirectionalPathIntegrator *CreateBidirectionalPathSurfaceIntegrator(const ParamSet &params);

#endif