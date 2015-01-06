BidirectionalPathTracer_PBRT-v2
===============================

Add BDPT to PBRT-v2:

- Put the .cpp and .h file on the scr/integrators directory

- Add the files to the project

- Open file core/api.cpp and look for the MakeSurfaceIntegrator function

- Add this code to the MakeSurfaceIntegrator funcion: 
    else if (name == "bidirectionalpath") si = CreateBidirectionalPathSurfaceIntegrator(paramSet);

- Build the project

NOTE: This algorithm has some errors. Feel free to correct them if you can. And then tell me :P
