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
uniform float u_noise_freq; 
uniform float u_step_length; 

out vec4 FragColor;

float rand(vec2 n); 
float noise(vec3 p); 
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
	float threshold_exp = 100000000.0; 
	//float inv_absortion_coef = 1.0/u_absortion_coef; 
	float optical_thickness = 0; 
	highp int res = int(u_noise_freq); 
	res = 8; 


	while( (num_step + 0.5) * u_step_length < inner_dist && 
			u_absortion_coef * optical_thickness < threshold_exp ) {

		vec3 curren_position = original_pos + ray_dir * num_step * u_step_length; 

		float increase_opt = u_step_length * noise(curren_position); 
		optical_thickness += increase_opt; 

		num_step += 1.0; 
	}
	
	optical_thickness = optical_thickness * u_absortion_coef; 

    float transmitansse = exp(-optical_thickness); 

    vec4 ret = vec4(u_color.xyz, u_color.w * (1.0 - transmitansse)); 
	//vec3 rand_vec = vec3(noise(v_world_position.xyz), noise(v_world_position.yzx), noise(v_world_position.zxy)); 
    //rand_vec = normalize(rand_vec); 
	//ret = vec4(rand_vec, 1); 

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

float mod289(float x){return x - floor(x * (1.0 / 289.0)) * 289.0;}
vec4 mod289(vec4 x){return x - floor(x * (1.0 / 289.0)) * 289.0;}
vec4 perm(vec4 x){return mod289(((x * 34.0) + 1.0) * x);}

float noise(vec3 p){
    vec3 a = floor(p);
    vec3 d = p - a;
    d = d * d * (3.0 - 2.0 * d);

    vec4 b = a.xxyy + vec4(0.0, 1.0, 0.0, 1.0);
    vec4 k1 = perm(b.xyxy);
    vec4 k2 = perm(k1.xyxy + b.zzww);

    vec4 c = k2 + a.zzzz;
    vec4 k3 = perm(c);
    vec4 k4 = perm(c + 1.0);

    vec4 o1 = fract(k3 * (1.0 / 41.0));
    vec4 o2 = fract(k4 * (1.0 / 41.0));

    vec4 o3 = o2 * d.z + o1 * (1.0 - d.z);
    vec2 o4 = o3.yw * d.x + o3.xz * (1.0 - d.x);

    return o4.y * d.y + o4.x * (1.0 - d.y);
}

vec2 vec_3_to_2(vec3 v) {
	// return vec2(v.x + 0.33 * v.z, v.y - 0.7421 * v.z); 
	return vec2(v.x + cos(v.z), v.y + sin(v.z)); 
}


