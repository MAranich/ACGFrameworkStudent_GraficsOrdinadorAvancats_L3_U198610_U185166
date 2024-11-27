#version 450 core
//#version 130 core

#define PI 3.14159265358979323846

in vec3 v_normal;
in vec3 v_position;
in vec3 v_world_position;
in vec2 v_uv;
in vec4 v_color;
in vec4 gl_FragCoord; 
// Uniforms
uniform mat4 u_model;
uniform mat4 u_viewprojection;
uniform vec3 u_camera_pos; 
uniform vec3 u_local_camera_position; 
uniform vec4 u_color; 
uniform vec3 u_bg_color; 
uniform vec3 u_box_min; 
uniform vec3 u_box_max; 
uniform float u_threshold; 
uniform float u_step_length; 

uniform sampler3D u_texture; 



out vec4 FragColor;


vec2 intersectAABB(vec3 rayOrigin, vec3 rayDir, vec3 boxMin, vec3 boxMax); 
vec3 world_to_local(vec3 world, mat4 model); 

float get_density(vec3 curren_position_local); 
 

vec3 ray_dir; 

void main() {
    
    // Converting the coordinates to Local
    vec3 pos = world_to_local(v_world_position, u_model); 
    // Get the direction of the ray
    ray_dir = pos - u_local_camera_position; 
    ray_dir = normalize(ray_dir); 

    //Finding the intersections w the auxiliary geometry
    vec2 int_dist = intersectAABB(u_local_camera_position, ray_dir, u_box_min, u_box_max); 
    
    //Definition of t near and t far 
	float t_near = int_dist.x; 
	float t_far = int_dist.y; 
    //The distance inside
    float inner_dist = t_far - t_near; 

    if(inner_dist <= 0.0) {
        FragColor = vec4(u_bg_color, 1); 
        return; 
    }

    //Ray equation, integration
	//vec3 curren_position = u_local_camera_position + ray_dir * (t_near + u_step_length * (0.5 + num_step)); 
	vec3 original_pos = u_local_camera_position + ray_dir * (t_near + u_step_length * 0.5); 
    //Variables init
	float num_step = 0; 


	while((num_step + 0.5) * u_step_length < inner_dist) { //While we are inside

        //Ray equation
        vec3 curren_position = original_pos + ray_dir * num_step * u_step_length; 
        
        // get_density() is the cnoise or texture3D, or 1 (case homogeneous)
        // conversion to texture coord. 
        vec3 curren_position_texture = (curren_position + 1.0f) * 0.5f;
        float density = texture3D(u_texture, curren_position_texture).x; 

        if(u_threshold <= density) {
            FragColor = u_color; 
            return; 
        }

        num_step += 1.0; 
	}
	
    FragColor = vec4(u_bg_color, 1); 
}

// DEFINITION OF FUNCTIONS ---------------------------------------------------------------
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


/*
// Depending on the density we are using: Homogeneous, Heterogeneous (Noise), or Bunny
float get_density(vec3 curren_position_local) {
    float ret = 1.0f; //Homogeneous
    switch (u_source_density) {
        //case 0: ret = 1.0f; 
    case 1: //Heterogeneous (Noise)
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

*/
