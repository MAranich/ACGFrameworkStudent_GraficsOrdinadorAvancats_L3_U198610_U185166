#version 450 core
//#version 130 core

#define PI 3.14159265358979323846
#define HEMISPHERICAL_NORM_CONSTANT 1.0/(4.0 * PI)
#define MAX_LIGHT 8

in vec3 v_normal;
in vec3 v_position;
in vec3 v_world_position;
in vec2 v_uv;
in vec4 v_color;
// Uniforms
uniform mat4 u_model;
uniform mat4 u_viewprojection;
uniform vec3 u_camera_pos; 
uniform vec3 u_local_camera_position; 
uniform vec4 u_color; 
uniform vec3 u_bg_color; 
uniform vec3 u_box_min; 
uniform vec3 u_box_max; 
uniform float u_absortion_coef; 
uniform float u_scattering_coef; 
uniform float u_step_length; 
uniform float u_scale;
uniform float u_detail;
uniform int u_source_density; 
uniform sampler3D u_texture; 
//Light
uniform int u_num_lights; // invariant: 0 <= u_num_lights <= MAX_LIGHT
uniform float u_light_intensity[MAX_LIGHT]; 
uniform vec3 u_light_color[MAX_LIGHT]; 
//uniform vec3 u_light_pos[MAX_LIGHT]; 
uniform vec3 u_light_pos_local[MAX_LIGHT]; 

uniform int u_num_scattering_steps; 
uniform int u_phase_fn; 
uniform float u_g; 

/*
    u_source_density

    0: uniform
    1: random
    2: bunny
*/



out vec4 FragColor;


vec2 intersectAABB(vec3 rayOrigin, vec3 rayDir, vec3 boxMin, vec3 boxMax); 
vec3 world_to_local(vec3 world, mat4 model); 

float get_density(vec3 curren_position_local); 

vec3 get_in_scattering(vec3 current_pos); 
float phase_function(float theta); 


#define MAX_OCTAVES 16

float hash1( float n ); 
float noise( vec3 x ); 
float fractal_noise( vec3 P, float detail ); 
float cnoise( vec3 P, float scale, float detail ); 


vec3 ray_dir; 

void main() {

    // FragColor = vec4(u_light_pos_local[0], 1.0f); 

    vec3 pos = world_to_local(v_world_position, u_model); 
    ray_dir = pos - u_local_camera_position; 
    ray_dir = normalize(ray_dir); 

    vec2 int_dist = intersectAABB(u_local_camera_position, ray_dir, u_box_min, u_box_max); 

	float t_near = int_dist.x; 
	float t_far = int_dist.y; 
    float inner_dist = t_far - t_near; 

    if(inner_dist <= 0.0) {
        FragColor = vec4(0, 0, 0, 0); 
        return; 
    }

	//vec3 curren_position = u_local_camera_position + ray_dir * (t_near + u_step_length * (0.5 + num_step)); 
	vec3 original_pos = u_local_camera_position + ray_dir * (t_near + u_step_length * 0.5); //Ray equation, integration
    //Variables init
	float num_step = 0; 
	float optical_thickness = 0; 
    vec3 pixel_color = vec3(0, 0, 0); 


	while((num_step + 0.5) * u_step_length < inner_dist) { //while we are inside

        vec3 curren_position = original_pos + ray_dir * num_step * u_step_length; //Ray equation
        
        // get_density() is the cnoise or texture3D, or 1 (case homogeneous)
        float current_absortion = get_density(curren_position); 

        // density * (ut = ua+ us)
        float transformed_absortion = current_absortion * (u_absortion_coef + u_scattering_coef);  
 

        //Riemann sum for heterogeneous
        optical_thickness += u_step_length * transformed_absortion;  

        //Transimittance (T(t',t))
        float atenuation = exp(-optical_thickness); 

        // ABSORPTION + OUT-SCATTERING TERM (using Le = u_color.xyz - EMISSION)
        pixel_color += u_color.xyz * transformed_absortion * atenuation * u_step_length; 


        // IN-SCATTERING TERM
        vec3 scattering_color = get_in_scattering(curren_position);
        float transformed_absorption_scat = current_absortion * u_scattering_coef; //us*density
        pixel_color += scattering_color * transformed_absorption_scat * atenuation * u_step_length;  

        num_step += 1.0; 
	}
	
    // Use the optical thickness accumulated computed to compute the Transmittance
    float transmitansse = exp(-optical_thickness); 
    
    //ret: the final L(t)
    vec4 ret = vec4(u_bg_color * transmitansse + pixel_color, 1.0); 

    FragColor = ret; 
}


// adapted from intersectCube in https://github.com/evanw/webgl-path-tracing/blob/master/webgl-path-tracing.js
// compute the near and far intersections of the cube (stored in the x and y components) using the slab method
// no intersection means vec.y - vec.x < 0 
vec2 intersectAABB(vec3 rayOrigin, vec3 rayDir, vec3 boxMin, vec3 boxMax) {
    vec3 tMin = (boxMin - rayOrigin) / rayDir;
    vec3 tMax = (boxMax - rayOrigin) / rayDir;
    vec3 t1 = min(tMin, tMax);
    vec3 t2 = max(tMin, tMax);
    float tNear = max(max(t1.x, t1.y), t1.z);
    float tFar = min(min(t2.x, t2.y), t2.z);
    return vec2(tNear, tFar);
};

// auxiliar function to transform world coordiantes to local ones
vec3 world_to_local(vec3 world, mat4 model) {

    vec4 tmp = vec4(world, 1.0); 
    tmp = model * tmp; 
    float inv_w = 1.0 / tmp.w; 

    return inv_w * tmp.xyz; 
}

vec3 get_in_scattering(vec3 origin_pos) { // we have to obtain the light from the main ray and the shadow ray 
    // origin_pos is local coords

    vec3 ret = vec3(0, 0, 0); 
    float inv_num_scattering_steps = 1.0f / float(u_num_scattering_steps); 

    for(int l = 0; l < u_num_lights; l++) {

        vec3 shadow_ray = normalize(u_light_pos_local[l] - origin_pos); 
        vec2 intersect = intersectAABB(origin_pos, shadow_ray, u_box_min, u_box_max); 
        //float final = max(intersect.x, intersect.y); 
        float final = intersect.y;

        //we divide final.y by the number of steps inside
        float step_length = final * inv_num_scattering_steps;  

        float optical_thickness = 0.0f; 

        for(int i = 0; i < u_num_scattering_steps; i++) {

            vec3 current_pos = origin_pos + shadow_ray * step_length * i; 
            float transformed_absortion = get_density(current_pos) * u_absortion_coef;//density*ua
            optical_thickness += transformed_absortion * step_length; // Riemann sum
        }

        float transmittance = exp(-optical_thickness); 
        //phase_function
        vec3 light_irradiance = u_light_color[l] * u_light_intensity[l]; 
        float angle = acos(dot(ray_dir, shadow_ray)); // we obtain the angle between the main ray and the shadow ray
        vec3 in_scattered_light = phase_function(angle) * light_irradiance; // we obtained the in-scattered light Ls

        vec3 current_pixel_col = in_scattered_light * transmittance; // we multiply Ls*T (transmittace)
        ret += current_pixel_col; 
    }

    return ret; 
}



float get_density(vec3 curren_position_local) {
    float ret = 1.0f; 
    switch (u_source_density) {
        //case 0: ret = 1.0f; 
    case 1:
        //random
        ret = cnoise(curren_position_local, u_scale, u_detail); 
        break; 
    case 2:
         //bunny

         // tex -> local : uv * 2 - 1
         // local -> tex : loc + 1 / 2

         vec3 curren_position_texture = (curren_position_local + 1.0f) * 0.5f;
         ret = texture3D(u_texture, curren_position_texture).x; 
         break; 

    default: 
        // also case = 0
        // ret = 1.0f; 
        break; 
        
    }

    return ret; 
}

float phase_function(float theta) { 
    float g = u_g; 
    float ret = 0.0f; 
    switch (u_phase_fn) {
        //case 0: //isotropic
    case 1: 
        float cos_theta = cos(theta); 
        float denominator = 1.0f + g * g - 2.0f * g * cos_theta; 
        float den_3_2 = sqrt(denominator * denominator * denominator); 
        ret = HEMISPHERICAL_NORM_CONSTANT * (1.0f - g * g) / den_3_2; 
        break; 
    default : 
        //case 0 = isotropic
        ret = HEMISPHERICAL_NORM_CONSTANT; 
    
    }

    return ret; 
}


// NOISE FUNCTIONS ///////////////////////////////////////////////////////

float hash1( float n )
{
    return fract( n*17.0*fract( n*0.3183099 ) );
}

float noise( vec3 x )
{
    vec3 p = floor(x);
    vec3 w = fract(x);
    
    vec3 u = w*w*w*(w*(w*6.0-15.0)+10.0);
    
    float n = p.x + 317.0*p.y + 157.0*p.z;
    
    float a = hash1(n+0.0);
    float b = hash1(n+1.0);
    float c = hash1(n+317.0);
    float d = hash1(n+318.0);
    float e = hash1(n+157.0);
    float f = hash1(n+158.0);
    float g = hash1(n+474.0);
    float h = hash1(n+475.0);

    float k0 =   a;
    float k1 =   b - a;
    float k2 =   c - a;
    float k3 =   e - a;
    float k4 =   a - b - c + d;
    float k5 =   a - c - e + g;
    float k6 =   a - b - e + f;
    float k7 = - a + b + c - d + e - f - g + h;

    return -1.0+2.0*(k0 + k1*u.x + k2*u.y + k3*u.z + k4*u.x*u.y + k5*u.y*u.z + k6*u.z*u.x + k7*u.x*u.y*u.z);
}

float fractal_noise( vec3 P, float detail )
{
    float fscale = 1.0;
    float amp = 1.0;
    float sum = 0.0;
    float octaves = clamp(detail, 0.0, 16.0);
    int n = int(octaves);

    for (int i = 0; i <= MAX_OCTAVES; i++) {
        if (i > n) continue;
        float t = noise(fscale * P);
        sum += t * amp;
        amp *= 0.5;
        fscale *= 2.0;
    }

    return sum;
}

float cnoise( vec3 P, float scale, float detail )
{
    P *= scale;
    return clamp(fractal_noise(P, detail), 0.0, 1.0);
}
