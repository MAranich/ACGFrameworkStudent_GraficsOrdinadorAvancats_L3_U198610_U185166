#version 450 core
//#version 130 core

#define PI 3.14159265358979323846
#define INV_PI 1.0 / 3.14159265358979323846
#define h 0.0002
//#define h 0.0005
#define TO_RGB 1.0/255.0
#define MAX_LIGHT 8

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

uniform int u_source_density; 

/*
    u_source_density

    0: uniform
    1: random
    2: bunny
*/

uniform float u_scale;
uniform float u_detail;

uniform bool u_use_jittering; 


//////////

uniform int u_mode; 
//uniform vec3 u_diffuse_color; 
uniform vec3 u_specular_color; 
uniform int u_shininess; 

//Light
uniform int u_num_lights; // invariant: 0 <= u_num_lights <= MAX_LIGHT
uniform float u_light_intensity[MAX_LIGHT]; 
uniform vec3 u_light_color[MAX_LIGHT]; 
uniform vec3 u_light_pos_local[MAX_LIGHT]; 


out vec4 FragColor;


vec2 intersectAABB(vec3 rayOrigin, vec3 rayDir, vec3 boxMin, vec3 boxMax); 
vec3 world_to_local(vec3 world, mat4 model); 

float get_density(vec3 curren_position_local); 
 
vec3 get_radiance(vec3 curren_position); 
vec3 get_gradient(vec3 current_position); 
vec3 reflectance_phong(vec3 incident_dir, vec3 normal); 


#define MAX_OCTAVES 16

float hash1( float n ); 
float noise( vec3 x ); 
float fractal_noise( vec3 P, float detail ); 
float cnoise( vec3 P, float scale, float detail ); 
float rand(vec2 co); 


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

    //TASK 1: Using jittering
    if(u_use_jittering) {
        original_pos = u_local_camera_position + ray_dir * (t_near + u_step_length * rand(gl_FragCoord.xy)); //using rand() function
    } 
    
    //Variables init
	float num_step = 0; 


	while((num_step + 0.5) * u_step_length < inner_dist) { //While we are inside

        //Ray equation
        vec3 curren_position = original_pos + ray_dir * num_step * u_step_length; 
        
        // get_density() is the cnoise or texture3D, or 1 (case homogeneous)
        float density = get_density(curren_position); 
        // We cast the ray and march through the volume until we reach the threshold
        if(u_threshold <= density) {

            switch (u_mode) {
                // case 0 = deafult
                case 1: 
                    // normals
                    float minimum_normal_length = 0.01; 

                    vec3 normal = -get_gradient(curren_position); //Normal is the inverse of the gradient
                    float normal_len_sq = dot(normal, normal); 
                    if(normal_len_sq < minimum_normal_length) { 
                        normal = vec3(0); 
                    } else {
                        normal = normalize(normal); 
                    }
                    normal = (normal + 1.0) * 0.5; //Texture coordinates

                    FragColor = vec4(normal, 1); 
                    break; 

                case 2: 
                    // phong coloring
                    FragColor = vec4(get_radiance(curren_position), 1); 
                    break; 
                default : 
                    // isosurface
                    FragColor = u_color; 
                    break; 
            }

            /*

                if(mode == 0) {
                    // isosurface
                    FragColor = u_color; 
                } else if(mode == 1) {
                    //normal

                    float minimum_normal_length = 0.01; 

                    vec3 normal = -get_gradient(curren_position); //Normal is the inverse of the gradient
                    float normal_len_sq = dot(normal, normal); 
                    if(normal_len_sq < minimum_normal_length) { 
                        normal = vec3(0); 
                    } else {
                        normal = normalize(normal); 
                    }
                    normal = (normal + 1.0) * 0.5; //Texture coordinates

                    FragColor = vec4(normal, 1); 
                } else {
                    // phong coloring
                    FragColor = vec4(get_radiance(curren_position), 1); 
                } 
            */
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

vec3 get_radiance(vec3 curren_position) {
    
    // u_color:  is diffuse coefititient and ambient coefitient
    vec3 ret = vec3(0); 

    float minimum_normal_length = 0.01; 

    vec3 normal = -get_gradient(curren_position); 
    float normal_len_sq = dot(normal, normal); 
    if(normal_len_sq < minimum_normal_length) { 
        normal = vec3(0); 
    } else {
        normal = normalize(normal); 
    }

    for(int l = 0; l < u_num_lights; l++) { 
        
        vec3 incident_dir = u_light_pos_local[l] - curren_position; //wi
        if (dot(incident_dir, normal) <= 0.0f) { 
            // visibility term is 0
            continue; 
        }

        vec3 light = u_light_color[l] * u_light_intensity[l]; //Li 
        
        incident_dir = normalize(incident_dir); 
        float dot_prod = dot(incident_dir, normal);  //wi·N
        
        vec3 brdf = reflectance_phong(incident_dir, normal); //fr

        ret += light * brdf * dot_prod; //V*Li*fr*(wi·N)

    }

    float ambient = 0.1f; 
    ret += u_color.xyz * ambient;  // L(wo) = V*Li*fr*(wi·N) + ambient

    return ret; 
        
    

}

vec3 get_gradient(vec3 current_position_local) {

    vec3 vec_x = vec3(h, 0, 0); 
    vec3 vec_y = vec3(0, h, 0); 
    vec3 vec_z = vec3(0, 0, h); 

    // position in texture coordinates
    vec3 position = (current_position_local + 1) * 0.5; 
   
    float dx = texture3D(u_texture, position + vec_x).x - texture3D(u_texture, position - vec_x).x; 
    float dy = texture3D(u_texture, position + vec_y).x - texture3D(u_texture, position - vec_y).x; 
    float dz = texture3D(u_texture, position + vec_z).x - texture3D(u_texture, position - vec_z).x; 
    
    float normalizer = 1.0f / (2.0f * h); 
    return normalizer * vec3(dx, dy, dz);  //Gradient
}

vec3 reflectance_phong(vec3 incident_dir, vec3 normal) {
    // w_i = incident_dir
    // w_0 = ray_dir (direction from point to camera * +-1)
    // w_r = (reflection ray_dir)

    ////////// DIFFUSE ////////////////
    // u_color is the diffuse coefitient
    vec3 diffuse = u_color.xyz * INV_PI; 

    ////////// SPECULAR ////////////////

    vec3 reflected_dir = 2.0 * dot(incident_dir, normal) * normal - incident_dir;  
    float dot_prod = dot(reflected_dir, ray_dir); 
    vec3 specular = u_specular_color * (2.0f * PI / (1.0f + u_shininess)) * pow(dot_prod, u_shininess); 


    return diffuse + specular; 
}

// NOISE FUNCTIONS ///////////////////////////////////////////////////////

float rand(vec2 co) {
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

float hash1( float n ) {
    return fract( n * 17.0 * fract( n * 0.3183099 ) );
}

float noise( vec3 x ) {
    vec3 p = floor(x);
    vec3 w = fract(x);
    
    vec3 u = w * w * w * (w * (w * 6.0 - 15.0) + 10.0);
    
    float n = p.x + 317.0*p.y + 157.0*p.z;
    
    float a = hash1(n+0.0);
    float b = hash1(n+1.0);
    float c = hash1(n+317.0);
    float d = hash1(n+318.0);
    float e = hash1(n+157.0);
    float f = hash1(n+158.0);
    float g = hash1(n+474.0);
    float hh = hash1(n+475.0);

    float k0 =   a;
    float k1 =   b - a;
    float k2 =   c - a;
    float k3 =   e - a;
    float k4 =   a - b - c + d;
    float k5 =   a - c - e + g;
    float k6 =   a - b - e + f;
    float k7 = - a + b + c - d + e - f - g + hh;

    return -1.0+2.0*(k0 + k1*u.x + k2*u.y + k3*u.z + k4*u.x*u.y + k5*u.y*u.z + k6*u.z*u.x + k7*u.x*u.y*u.z);
}

float fractal_noise( vec3 P, float detail ) {
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

float cnoise( vec3 P, float scale, float detail ) {
    P *= scale;
    return clamp(fractal_noise(P, detail), 0.0, 1.0);
}
