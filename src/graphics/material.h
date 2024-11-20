#pragma once

#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/matrix.hpp>

#include "../framework/camera.h"
#include "mesh.h"
#include "texture.h"
#include "shader.h"

// \libraries\easyVDB\src
#include "openvdbReader.h"
#include "bbox.h"


enum VolumeDensityMode { Homogeneus, Noise, Bunny };
enum PhaseFunction { Isotropic, HenyeyGreenstein, Cardioid };

class Material {
public:

	Shader* shader = NULL;
	Texture* texture = NULL;
	glm::vec4 color;

	virtual void setUniforms(Camera* camera, glm::mat4 model) = 0;
	virtual void render(Mesh* mesh, glm::mat4 model, Camera* camera) = 0;
	virtual void renderInMenu() = 0;
};

class FlatMaterial : public Material {
public:

	FlatMaterial(glm::vec4 color = glm::vec4(1.f));
	~FlatMaterial();

	void setUniforms(Camera* camera, glm::mat4 model);
	void render(Mesh* mesh, glm::mat4 model, Camera* camera);
	void renderInMenu();
};

class WireframeMaterial : public FlatMaterial {
public:

	WireframeMaterial();
	~WireframeMaterial();

	void render(Mesh* mesh, glm::mat4 model, Camera* camera);
};

class StandardMaterial : public Material {
public:

	bool first_pass = false;

	bool show_normals = false;
	Shader* base_shader = NULL;
	Shader* normal_shader = NULL;

	StandardMaterial(glm::vec4 color = glm::vec4(1.f));
	~StandardMaterial();

	void setUniforms(Camera* camera, glm::mat4 model);
	void render(Mesh* mesh, glm::mat4 model, Camera* camera);
	void renderInMenu();


};


class VolumeMaterial : public Material {
	public:
		~VolumeMaterial();
		VolumeMaterial(glm::vec4 color);

		float absortion_coefitient;
		float scattering_coefitient;
		VolumeDensityMode density_mode; 
		float density_plus;

		std::vector<Shader*> shader_list;
		float step_length = 0.05f; 
		float scale = 2.209f;
		float detail = 5.0f;
		int num_scatter_steps = 5;

		PhaseFunction phase_function; 
		float g_coef; 

		void setUniforms(Camera* camera, glm::mat4 model);
		void render(Mesh* mesh, glm::mat4 model, Camera* camera);
		void renderInMenu();

		void loadVDB(std::string file_path);

		void estimate3DTexture(easyVDB::OpenVDBReader* vdbReader);

};