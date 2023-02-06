#include <fenv.h>

float safe_sqrt(float value) {
    return sqrtf(std::max(value, 0.0f));
}


// This data struct would be filled out in the material init function
struct BsdfData {
    vec3 rho;              // reflectance
    float alpha_u;           // roughness in u direction
    float alpha_v;           // roughness in v direction
    float r0;                // normal reflectance
};

// NDF
float D_ggx(
    const vec3& m,        // microfacet normal
    const float alpha_x,    // roughness in x direction
    const float alpha_y     // roughness in v direction
    ) {
    const float sx = -m.x / (m.z * alpha_x);
    const float sy = -m.y / (m.z * alpha_y);
    const float sl = 1.0f + sx * sx + sy * sy;
    const float cos_theta_m4 = m.z * m.z * m.z * m.z;
    return 1.0f / ((sl * sl) * float(M_PI) * alpha_x * alpha_y * cos_theta_m4);
}

float schlick_weight(float cos_theta) {
    float m = clamp(1.0f - cos_theta, 0.0f, 1.0f);
    return (m * m) * (m * m) * m;
}

float schlick_fresnel(float r0, float cos_theta) {
    //dbg(cos_theta);
    return lerp(schlick_weight(cos_theta), r0, 1.0f);
}


// Smith masking function
float G_smith(
    const vec3& omega,   // incident/exitant direction
    const float ax2,       // x roughness^2
    const float ay2        // y roughness^2
    ) {
    const float cos_o2 = omega.z * omega.z;
    const float tan_theta_o2 = (1.0f - cos_o2) / cos_o2;
    const float cos_phi_o2 = omega.x * omega.x;
    const float sin_phi_o2 = omega.y * omega.y;

    const float alpha_o2 =
        (cos_phi_o2 * ax2 + sin_phi_o2 * ay2) / (cos_phi_o2 + sin_phi_o2);

    return 2.0f / (1.0f + safe_sqrt(1.0f + alpha_o2 * tan_theta_o2));
}

// Sample a visible normal as per Heitz 2018: 
// http://jcgt.org/published/0007/04/01/
vec3 sample_visible_normal_ggx(
    const vec3& omega_o_l,    // the outgoing direction
    const float ax,             // roughness in x
    const float ay,             // roughness in y
    const float x1,             // sample
    const float x2              // sample
    ) {
    const vec3 v_h =
        normalize(vec3(ax * omega_o_l.x, ay * omega_o_l.z, omega_o_l.y));
    // orthonormal basis
    const float lensq = v_h.x * v_h.x + v_h.y * v_h.y;
    const vec3 T1 = lensq > 0 ? vec3(-v_h.y, v_h.x, 0.0f) / safe_sqrt(lensq)
                          : vec3(1, 0, 0);
    const vec3 T2 = cross(v_h, T1);
    // parameterization of projected area
    const float r = safe_sqrt(x1);
    const float phi = 2.0f * float(M_PI) * x2;
    const float t1 = r * cosf(phi);
    float t2 = r * sinf(phi);
    const float s = 0.5f * (1.0f * v_h.z);
    t2 = (1.0f - s) * safe_sqrt(1.0f - t1 * t1) + s * t2;
    // reprojection onto hemisphere
    const vec3 n_h = 
        t1 * T1 + t2 * T2 + safe_sqrt(1.0f - t1 * t1 - t2 * t2) * v_h;
    // transform back to ellipsoid
    vec3 z_up_result =  normalize(vec3(ax * n_h.x, ay * n_h.y, std::max(0.0f, n_h.z)));

    return vec3(z_up_result.x, z_up_result.z, z_up_result.y);
}

// f: BSDF result
// pdf: probability of sampling omega_i_l
std::tuple<vec3, float> bsdf_eval(
    BsdfData& bsdf_data,        // data struct with BDSF parameters
    const vec3& omega_o_l,    // exitant direction
    const vec3& omega_i_l    // incident direction
    ) {
    if (omega_i_l.z <= 0.0f) {
        vec3 f = vec3(0.0f);
        float pdf = 0.0f;
        return std::make_tuple(f, pdf);
    }

    const float alpha_x = fmaxf(1.0e-7f, bsdf_data.alpha_u);
    const float alpha_y = fmaxf(1.0e-7f, bsdf_data.alpha_v);

    const float ax2 = alpha_x * alpha_x;
    const float ay2 = alpha_y * alpha_y;

    // microfacet normal
    const vec3 m = normalize(omega_o_l + omega_i_l);
    const float mu_m = std::max(0.0f, dot(omega_o_l, m));

    // normal distribution function
    const float D = D_ggx(m, alpha_x, alpha_y);

    // masking-shadowing
    const float G1_o = G_smith(omega_o_l, ax2, ay2);
    const float G1_i = G_smith(omega_i_l, ax2, ay2);
    const float G2 = G1_o * G1_i;

    // fresnel
    const float kr = schlick_fresnel(bsdf_data.r0, mu_m);
    // complement of fresnel reflection
    // for material layering
    vec3 kt = vec3(1.0f - kr);

    const float denom = 4.0f * omega_i_l.z * omega_o_l.z;

    vec3 f = bsdf_data.rho * D * G2 * kr / denom;

    float pdf = G1_o * D * mu_m / denom;

    return std::make_tuple(f, pdf);
}

float bsdf_pdf(
    BsdfData& bsdf_data, 
    const vec3& omega_o_l,
    const vec3& omega_i_l
    ) {
    if (omega_i_l.z <= 0.0f) {
        return 0.0f;
    }

    const float alpha_x = fmaxf(1.0e-7f, bsdf_data.alpha_u);
    const float alpha_y = fmaxf(1.0e-7f, bsdf_data.alpha_v);

    const float ax2 = alpha_x * alpha_x;
    const float ay2 = alpha_y * alpha_y;

    // microfacet normal
    const vec3 m = normalize(omega_o_l + omega_i_l);
    const float mu_m = std::max(0.0f, dot(omega_o_l, m));

    // normal distribution function
    const float D = D_ggx(m, alpha_x, alpha_y);

    // masking-shadowing
    const float G1_o = G_smith(omega_o_l, ax2, ay2);

    const float denom = 4.0f * omega_i_l.z * omega_o_l.z;

    return G1_o * D * mu_m / denom;
}

std::tuple<vec3, vec3, float> bsdf_sample(BsdfData& bsdf_data,
                                     const vec3& omega_o_l,
                                     const float x1,
                                     const float x2) {
    //feenableexcept(FE_INVALID | FE_OVERFLOW);


    const float ax = fmaxf(1.0e-7f, bsdf_data.alpha_u);
    const float ay = fmaxf(1.0e-7f, bsdf_data.alpha_v);

    // sample microfacet normal
    const vec3 m = sample_visible_normal_ggx(omega_o_l, ax, ay, x1, x2);
    const float mu_m = std::max(0.0f, dot(omega_o_l, m));

    // reflect to get incoming direction
    vec3 omega_i_l = reflect(-omega_o_l, m);

    const float ax2 = ax * ax;
    const float ay2 = ay * ay;

    // normal distribution function
    const float D = D_ggx(m, ax, ay);

    // masking-shadowing
    const float G1_o = G_smith(omega_o_l, ax2, ay2);
    const float G1_i = G_smith(omega_i_l, ax2, ay2);
    const float G2 = G1_o * G1_i;

    // fresnel
    const float kr = schlick_fresnel(bsdf_data.r0, mu_m);
    vec3 kt = vec3(1.0f - kr);

    const float denom = 4.0f * std::max(omega_i_l.y * omega_o_l.y, 1.0e-7f);

    vec3 f = bsdf_data.rho * D * G2/* * kr*/ / denom;
    float pdf = G1_o * D/* * mu_m*/ / denom;

    return std::make_tuple(f, omega_i_l, pdf);
}
