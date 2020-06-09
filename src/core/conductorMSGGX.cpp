
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// core/reflection.cpp*

#include "reflection.h"
#include "spectrum.h"
#include "sampler.h"
#include "sampling.h"
#include "interpolation.h"
#include "scene.h"
#include "interaction.h"
#include "stats.h"
#include "roughGGX.h"
#include <stdarg.h>

namespace pbrt {


Spectrum ConductorMSReflection::f(const Vector3f &wo, const Vector3f &wi) const {
    if (!SameHemisphere(wo, wi)) return Spectrum(0.f);
    Float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
    // Handle degenerate cases for microfacet reflection
    if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0.);
    //Vector3f wh = wi + wo;
    //if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.);


    Spectrum res = 2.0f * eval_conductor(wo, wi, alpha_x, alpha_y, m_eta, m_k, m_scatteringOrderMax)/cosThetaI;
    //Spectrum res = 2.0f * eval_conductor(wo, wi, alpha_x, alpha_y, m_eta, m_k, m_scatteringOrderMax);
    /*
    Spectrum res = (generateRandomNumber() > 0.5f) ?
                    2.0f * eval_conductor(wo, wi, alpha_x, alpha_y, m_eta, m_k, m_scatteringOrderMax)/cosThetaI:
                    2.0f * eval_conductor(wi, wo, alpha_x, alpha_y, m_eta, m_k, m_scatteringOrderMax)/cosThetaO;
    */
    return res; 

}

std::string ConductorMSReflection::ToString() const {
    return std::string("[ConductorMSReflection m_eta: ") + m_eta.ToString() +
           std::string(" m_k: ") + m_k.ToString() +
           StringPrintf(" alpha_x: %f; alpha_y: %f ]", alpha_x, alpha_y);
} 


Spectrum ConductorMSReflection::Sample_f(const Vector3f &wo, Vector3f *wi,
                                        const Point2f &u, Float *pdf,
                                        BxDFType *sampledType) const {

/*
    bool flip = wo.z < 0;
    Vector3f wh = TrowbridgeReitzSample(flip ? -wo : wo, alpha_x, alpha_y, u[0], u[1]);
    if (flip) wh = -wh;
    if (Dot(wo, wh) < 0) return 0.;   // Should be rare
    *wi = Reflect(wo, wh);
    if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);
    Spectrum energy = f(wo, *wi);
    if (pdf) *pdf = Pdf(wo, *wi);
    return energy;	
 
*/ 
    if (!useIS) {
        *wi = CosineSampleHemisphere(u);
        if (wo.z < 0) wi->z *= -1;
        *pdf = Pdf(wo, *wi);
        return f(wo, *wi);
    }

    Spectrum energy(1.0f); 
    Vector3f wi_sampled = sample_conductor(wo, alpha_x, alpha_y, m_eta, m_k, m_scatteringOrderMax, energy);    
    *wi = wi_sampled;
    if (!SameHemisphere(wo, *wi) || wi_sampled.z == 0) return Spectrum(0.f);
    Float thePdf = Pdf(wo, *wi);
    if (pdf) *pdf = thePdf;
    //energy *= (thePdf/fabs(wi_sampled.z));
    return energy  * thePdf / AbsCosTheta(wi_sampled);

}

Float ConductorMSReflection::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    if (!useIS) return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;

    if (!SameHemisphere(wo, wi)) return 0.0f;

    Vector3f wh = Normalize(wo+wi);

    RayInfo ray;
    ray.updateDirection(wo, alpha_x, alpha_y);

    /*
    if (sampleVisibleArea)
        return D(wh) * G1(wo) * AbsDot(wo, wh) / AbsCosTheta(wo);
    else
        return D(wh) * AbsCosTheta(wh);
    distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));
    */


                // single-scattering PDF + diffuse
                // otherwise too many fireflies due to lack of multiple-scattering PDF
                // (MIS works even if the PDF is wrong and not normalized)
    //float scale = InvPi;
    float scale = diffPdfScale;
    float diffusePdf = AbsCosTheta(wi) * scale;
    //return D_ggx(wh, alpha_x, alpha_y) / (1.0f + ray.Lambda) / (4.0f * AbsCosTheta(wo)) + diffusePdf;

    float p = distribution->Pdf(wo, wh) / (4 * Dot(wo, wh))  + diffusePdf;
    return p;
}

}  // namespace pbrt
