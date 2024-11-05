#version 450 core
//#version 130 core

#define PI 3.14159265358979323846

in vec3 v_normal;
in vec3 v_position;
in vec3 v_world_position;
in vec2 v_uv;
in vec4 v_color;

uniform mat4 u_model;
uniform mat4 u_viewprojection;
uniform vec3 u_camera_pos; 
uniform vec3 u_local_camera_position; 
uniform vec4 u_color; 
uniform vec3 u_box_min; 
uniform vec3 u_box_max; 
uniform float u_absortion_coef; 
uniform float u_step_length; 
uniform float u_scale;
uniform float u_detail;
uniform vec3 u_bg_color;

out vec4 FragColor;



float rand(vec2 n); 
float noise(vec3 p); 
vec2 vec_3_to_2(vec3 x); 
vec2 intersectAABB(vec3 rayOrigin, vec3 rayDir, vec3 boxMin, vec3 boxMax); 

vec3 world_to_local(vec3 world, mat4 model); 

// NOISE FUNCTIONS ///////////////////////////////////////////////////////
#define MAX_OCTAVES 16

float hash1( float n ); 
float noise( vec3 x ); 
float fractal_noise( vec3 P, float detail ); 
float cnoise( vec3 P, float scale, float detail ); 

//////////////////////////////////////////////////////////////////////////////////////////

void main() {
    
    vec3 pos = world_to_local(v_world_position, u_model); 
    vec3 ray_dir = pos - u_local_camera_position; 
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
	vec3 original_pos = u_local_camera_position + ray_dir * (t_near + u_step_length * 0.5); 
	float num_step = 0; 
	float threshold_exp = 10.0; 
	float optical_thickness = 0; 


	while( (num_step + 0.5) * u_step_length < inner_dist && 
			u_absortion_coef * optical_thickness < threshold_exp ) {

		vec3 curren_position = original_pos + ray_dir * num_step * u_step_length; 

		float increase_opt = u_step_length * cnoise(curren_position,u_scale, u_detail);
		optical_thickness += increase_opt; 

		num_step += 1.0; 
	}
	

    optical_thickness = (optical_thickness * u_absortion_coef); 
    float transmittance = exp(-optical_thickness);

    vec4 ret = vec4(u_bg_color.xyz * transmittance, 1.0); 

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