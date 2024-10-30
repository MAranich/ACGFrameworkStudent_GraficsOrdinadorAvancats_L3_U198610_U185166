//#version 450 core
#version 130 core

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
uniform vec4 u_bg_color;

out vec4 FragColor;

///////////////// EXPERIMETAL
in vec4 gl_FragCoord; 

////////////////

vec2 intersectAABB(vec3 rayOrigin, vec3 rayDir, vec3 boxMin, vec3 boxMax); 

vec3 world_to_local(vec3 world, mat4 model); 

void main() {
    // u_camera_pos is in world space
    // v_position is in ???
    
    // Color is (0, 0, 0, 1) ???

    // vec3 ray_dir = v_position - u_camera_pos; 
    vec3 pos = world_to_local(v_world_position, u_model); 
    vec3 ray_dir = pos - u_local_camera_position; 
    ray_dir = normalize(ray_dir); 

    vec2 int_dist = intersectAABB(u_local_camera_position, ray_dir, u_box_min, u_box_max); 

    float inner_dist = int_dist.y - int_dist.x; 
    // ^positive if collision

    if(inner_dist <= 0.0) {
        FragColor = u_bg_color; 
        return; 
    }

    // inner_dist = 0 => exp() = 1
    // inner_dist = 1000 => exp() = 0
    float optical_thickness = exp(-inner_dist * u_absortion_coef); 

    vec4 ret = vec4(u_color.xyz, u_color.w * (1.0-optical_thickness)); 
    //ret = vec4(vec3(inner_dist * u_absortion_coef), 1); 
    //ret = u_color; 
    //ret = vec4(u_color.xyz, 0.5); 
    //ret = vec4(0, 0, 0, u_color.w * (1.0-optical_thickness)); 

    FragColor = ret; 
    
}


// adapted from intersectCube in https://github.com/evanw/webgl-path-tracing/blob/master/webgl-path-tracing.js
// compute the near and far intersections of the cube (stored in the x and y components) using the slab method
// no intersection means vec.x > vec.y (really tNear > tFar)
// no intersection means vec.y < vec.x 
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
    float inv_w = 1/tmp.w; 

    return inv_w * tmp.xyz; 
}

