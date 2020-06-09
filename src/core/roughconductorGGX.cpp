
#include "reflection.h"
#include "roughGGX.h"

namespace pbrt {

Vector3f samplePhaseFunction_conductor(const Vector3f& wi, const float alpha_x, const float alpha_y, const Spectrum& m_eta, const Spectrum& m_k, Spectrum& weight)
{
	const float U1 = generateRandomNumber();
	const float U2 = generateRandomNumber();

	// sample D_wi
	// stretch to match configuration with alpha=1.0	
	const Vector3f wi_11 = Normalize(Vector3f(alpha_x * wi.x, alpha_y * wi.y, wi.z));

	// sample visible slope with alpha=1.0
	Vector3f slope_11 = sampleP22_11(acosf(wi_11.z), U1, U2, alpha_x, alpha_y);

	// align with view direction
	const float phi = atan2(wi_11.y, wi_11.x);
	Vector3f slope(cosf(phi)*slope_11.x - sinf(phi)*slope_11.y, sinf(phi)*slope_11.x + cos(phi)*slope_11.y, 0);

	// stretch back
	slope.x *= alpha_x;
	slope.y *= alpha_y;

	// compute normal
	Vector3f wm;
	// if numerical instability
	if( (slope.x != slope.x) || !IsFiniteNumber(slope.x) ) 
	{
		if(wi.z > 0) wm = Vector3f(0.0f,0.0f,1.0f);
		else wm = Normalize(Vector3f(wi.x, wi.y, 0.0f));
	}
	else
		wm = Normalize(Vector3f(-slope.x, -slope.y, 1.0f));

	// reflect
	const Vector3f wo = -wi + 2.0f * wm * Dot(wi, wm);
	weight = FrConductor(Dot(wi, wm), m_eta, m_k);

	return wo;
}

Spectrum evalPhaseFunction_conductor(const RayInfo& ray, const Vector3f& wo, const float alpha_x, const float alpha_y, const Spectrum& m_eta, const Spectrum& m_k)
{
	if(ray.w.z > 0.9999f)
		return Spectrum(0.0f);

	// half vector 
	const Vector3f wh = Normalize(-ray.w+wo);
	if(wh.z < 0.0f)
		return Spectrum(0.0f);
	
	// projected area
	float projectedArea;
	if(ray.w.z < -0.9999f)
		projectedArea = 1.0f;
	else 
		projectedArea = ray.Lambda * ray.w.z;

	// value
	const Spectrum value = FrConductor(Dot(-ray.w, wh), m_eta, m_k) * std::max(0.0f, Dot(-ray.w, wh)) * D_ggx(wh, alpha_x, alpha_y) / 4.0f / projectedArea / Dot(-ray.w, wh);
	return value;
}

// MIS weights for bidirectional path tracing on the microsurface
float MISweight_conductor(const Vector3f& wi, const Vector3f& wo, const float alpha_x, const float alpha_y) 
{
	if(wi.x == -wo.x && wi.y == -wo.y && wi.z == -wo.z)
		return 1.0f;
	const Vector3f wh = Normalize(wi+wo);
	const float value = D_ggx( (wh.z>0) ? wh : -wh , alpha_x, alpha_y);
	return value;
}

Spectrum eval_conductor(const Vector3f& wi, const Vector3f& wo, const float alpha_x, const float alpha_y, const Spectrum& m_eta, const Spectrum& m_k, const int scatteringOrderMax)
{
	if(wi.z <= 0 || wo.z <= 0)
		return Spectrum(0.0f);

	// init
	RayInfo ray;
	ray.updateDirection(-wi, alpha_x, alpha_y);	
	ray.updateHeight(1.0f);
	Spectrum energy(1.0f);

	RayInfo ray_shadowing;
	ray_shadowing.updateDirection(wo, alpha_x, alpha_y);

	// eval single scattering	
	// half-vector
	const Vector3f wh = Normalize(wi+wo);
	const float D = D_ggx(wh, alpha_x, alpha_y);
	const float G2 = 1.0f / (1.0f + (-ray.Lambda-1.0f) + ray_shadowing.Lambda);
	Spectrum singleScattering = FrConductor(Dot(-ray.w, wh), m_eta, m_k)  *  D * G2 / (4.0f * wi.z);
    if (scatteringOrderMax == 1)
	  return 0.5 * singleScattering;
	
	// MIS weight 
	float wi_MISweight;

	// multiple scattering
	Spectrum multipleScattering(0.0f);
	
	// random walk
	int current_scatteringOrder = 0;	
	while(current_scatteringOrder < scatteringOrderMax)
	{
		// next height
		float U = generateRandomNumber();
		ray.updateHeight( sampleHeight(ray, U) );		
				
		// leave the microsurface?
		if( ray.h == std::numeric_limits<Float>::max() )
			break;
		else
			current_scatteringOrder++;

		// next event estimation 
		if( current_scatteringOrder > 1) // single scattering is already computed
		{
			Spectrum phasefunction = evalPhaseFunction_conductor(ray, wo, alpha_x, alpha_y, m_eta, m_k); 
			ray_shadowing.updateHeight(ray.h);
			float shadowing = ray_shadowing.G1;
			Spectrum I = energy * phasefunction * shadowing;

			// MIS
			const float MIS = wi_MISweight / ( wi_MISweight + MISweight_conductor(-ray.w, wo, alpha_x, alpha_y) );


			if ( IsFiniteNumber(I[0]) )
				multipleScattering += I * MIS;
		}

		// next direction
		Spectrum weight;
		ray.updateDirection(samplePhaseFunction_conductor(-ray.w, alpha_x, alpha_y, m_eta, m_k, weight), alpha_x, alpha_y);
		energy = energy * weight;
		ray.updateHeight(ray.h);

		if(current_scatteringOrder == 1)
			wi_MISweight = MISweight_conductor(wi, ray.w, alpha_x, alpha_y);

		// if NaN (should not happen, just in case)
		if( (ray.h != ray.h) || (ray.w.x != ray.w.x)) 
			return Spectrum(0.0f);
	}

	// 0.5f = MIS weight of singleScattering
	// multipleScattering already weighted by MIS
	return 0.5f*singleScattering + multipleScattering;
}

Vector3f sample_conductor(const Vector3f& wi, const float alpha_x, const float alpha_y, const Spectrum& m_eta, const Spectrum& m_k, const int scatteringOrderMax, Spectrum& energy)
{
	energy = Spectrum(1.0f);

	// init
	RayInfo ray;
	ray.updateDirection(-wi, alpha_x, alpha_y);
	ray.updateHeight(1.0f);
		
	// random walk
	int current_scatteringOrder = 0;
	while(true)
	{
		// next height
		float U = generateRandomNumber();
		ray.updateHeight( sampleHeight(ray, U) );		

		// leave the microsurface?
		if( ray.h == std::numeric_limits<Float>::max() )
			break;
		else
			current_scatteringOrder++;

		// next direction
		Spectrum weight;
		ray.updateDirection(samplePhaseFunction_conductor(-ray.w, alpha_x, alpha_y, m_eta, m_k, weight), alpha_x, alpha_y);
		energy = energy * weight;
		ray.updateHeight(ray.h);

		// if NaN (should not happen, just in case)
		if( (ray.h != ray.h) || (ray.w.x != ray.w.x)) 
		{
			energy = Spectrum(0.0f);
			return Vector3f(0,0,1);
		}

		if( current_scatteringOrder > scatteringOrderMax )
		{
			energy = Spectrum(0.0f);
			return Vector3f(0,0,1);
		}
	}

	return ray.w;
}


}//end namespace

