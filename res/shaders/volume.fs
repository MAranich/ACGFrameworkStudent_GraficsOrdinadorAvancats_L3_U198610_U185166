#version 450 core

in vec3 v_normal;
in vec3 v_position;
in vec3 v_world_position;
in vec2 v_uv;
in vec4 v_color;


uniform vec3 u_camera_pos; 
uniform vec3 u_box_min; 
uniform vec3 u_box_max; 
uniform float absortion_coefitient; 

out vec4 FragColor;

///////////////// EXPERIMETAL
in vec4 gl_FragCoord; 

////////////////

vec2 intersectAABB(vec3 rayOrigin, vec3 rayDir, vec3 boxMin, vec3 boxMax); 


void main() {
    
    vec3 ray_dir = v_position - u_camera_pos; 
    ray_dir = normalize(ray_dir); 

    vec2 int_dist = intersectAABB(u_camera_pos, ray_dir, u_box_min, u_box_max); 

    float inner_dist = int_dist.x - int_dist.y; 

    if(inner_dist <= 0.0) {
        FragColor = v_color; 
    }

    // vec4 bg_color = vec4(0.1, 0.1, 0.1, 1.0); 

    float optical_thickness = exp(-inner_dist * absortion_coefitient); 

    vec4 ret = vec4(v_color.xyz, v_color.w * (1.0 - optical_thickness)); 

    FragColor = ret; 
    
}


// adapted from intersectCube in https://github.com/evanw/webgl-path-tracing/blob/master/webgl-path-tracing.js

// compute the near and far intersections of the cube (stored in the x and y components) using the slab method
// no intersection means vec.x > vec.y (really tNear > tFar)
vec2 intersectAABB(vec3 rayOrigin, vec3 rayDir, vec3 boxMin, vec3 boxMax) {
    vec3 tMin = (boxMin - rayOrigin) / rayDir;
    vec3 tMax = (boxMax - rayOrigin) / rayDir;
    vec3 t1 = min(tMin, tMax);
    vec3 t2 = max(tMin, tMax);
    float tNear = max(max(t1.x, t1.y), t1.z);
    float tFar = min(min(t2.x, t2.y), t2.z);
    return vec2(tNear, tFar);
};
