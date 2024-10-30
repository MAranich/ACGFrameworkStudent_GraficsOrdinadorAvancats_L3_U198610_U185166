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
uniform float u_absortion_coef_mult; 
uniform float u_noise_freq; 
uniform float u_step_length; 

out vec4 FragColor;

float rand(vec2 n); 
float noise(vec2 p, float freq); 
float pNoise(vec2 p, int res); 
vec2 vec_3_to_2(vec3 x); 
vec2 intersectAABB(vec3 rayOrigin, vec3 rayDir, vec3 boxMin, vec3 boxMax); 

vec3 world_to_local(vec3 world, mat4 model); 

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

	vec3 original_pos = u_local_camera_position + (t_near + u_step_length * 0.5) * ray_dir; 
	float num_step = 0; 
	float threshold_exp = 10.0; 
	float optical_thickness = 0; 

	while( (num_step + 0.5) * u_step_length < inner_dist && optical_thickness < threshold_exp ) {

		vec3 curren_position = original_pos + ray_dir * num_step * u_step_length; 

		highp int res = int(u_noise_freq); 
		float increase_opt = u_step_length * pNoise(vec_3_to_2(curren_position), res); 
		optical_thickness += increase_opt; 

		num_step += 1.0; 
	}
	
	optical_thickness = optical_thickness * u_absortion_coef_mult; 

    float transmitansse = exp(-optical_thickness); 

    vec4 ret = vec4(u_color.xyz, u_color.w * (1.0-optical_thickness)); 

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


float rand(vec2 n) { 
	return fract(sin(dot(n, vec2(12.9898, 4.1414))) * 43758.5453);
}

float noise(vec2 p, float freq ){
	float unit = 1024/freq;
	//float unit = screenWidth/freq;
	vec2 ij = floor(p/unit);
	vec2 xy = mod(p,unit)/unit;
	//xy = 3.*xy*xy-2.*xy*xy*xy;
	xy = .5*(1.-cos(PI*xy));
	float a = rand((ij+vec2(0.,0.)));
	float b = rand((ij+vec2(1.,0.)));
	float c = rand((ij+vec2(0.,1.)));
	float d = rand((ij+vec2(1.,1.)));
	float x1 = mix(a, b, xy.x);
	float x2 = mix(c, d, xy.x);
	return mix(x1, x2, xy.y);
}

// Perlin noise
float pNoise(vec2 p, int res){
	float persistance = .5;
	float n = 0.;
	float normK = 0.;
	float f = 4.;
	float amp = 1.;
	int iCount = 0;
	for (int i = 0; i<50; i++){
		n += amp * noise(p, f);
		f*=2.;
		normK+=amp;
		amp*=persistance;
		if (iCount == res) break;
		iCount++;
	}
	float nf = n/normK;
	return nf*nf*nf*nf;
}

vec2 vec_3_to_2(vec3 v) {
	// return vec2(v.x + 0.33 * v.z, v.y - 0.7421 * v.z); 
	return vec2(v.x + cos(v.z), v.y + sin(v.z)); 
}
