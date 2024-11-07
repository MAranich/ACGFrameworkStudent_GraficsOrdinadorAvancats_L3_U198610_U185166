#include "material.h"

#include "application.h"

#include <istream>
#include <fstream>
#include <algorithm>
#define MAX_LIGHT 8

FlatMaterial::FlatMaterial(glm::vec4 color)
{
	this->color = color;
	this->shader = Shader::Get("res/shaders/basic.vs", "res/shaders/flat.fs");
}

FlatMaterial::~FlatMaterial() { }

void FlatMaterial::setUniforms(Camera* camera, glm::mat4 model)
{
	//upload node uniforms
	this->shader->setUniform("u_viewprojection", camera->viewprojection_matrix);
	this->shader->setUniform("u_camera_position", camera->eye);
	this->shader->setUniform("u_model", model);

	this->shader->setUniform("u_color", this->color);

	this->shader->setUniform("u_source_density", 0);
}

void FlatMaterial::render(Mesh* mesh, glm::mat4 model, Camera* camera)
{
	if (mesh && this->shader) {
		// enable shader
		this->shader->enable();

		// upload uniforms
		setUniforms(camera, model);

		// do the draw call
		mesh->render(GL_TRIANGLES);

		this->shader->disable();
	}
}

void FlatMaterial::renderInMenu()
{
	ImGui::ColorEdit3("Color", (float*)&this->color);
}

WireframeMaterial::WireframeMaterial()
{
	this->color = glm::vec4(1.f);
	this->shader = Shader::Get("res/shaders/basic.vs", "res/shaders/flat.fs");
}

WireframeMaterial::~WireframeMaterial() { }

void WireframeMaterial::render(Mesh* mesh, glm::mat4 model, Camera* camera)
{
	if (this->shader && mesh)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDisable(GL_CULL_FACE);

		//enable shader
		this->shader->enable();

		//upload material specific uniforms
		setUniforms(camera, model);

		//do the draw call
		mesh->render(GL_TRIANGLES);

		glEnable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
}

StandardMaterial::StandardMaterial(glm::vec4 color)
{
	this->color = color;
	this->base_shader = Shader::Get("res/shaders/basic.vs", "res/shaders/basic.fs");
	this->normal_shader = Shader::Get("res/shaders/basic.vs", "res/shaders/normal.fs");
	this->shader = this->base_shader;
}

StandardMaterial::~StandardMaterial() { }

void StandardMaterial::setUniforms(Camera* camera, glm::mat4 model)
{
	//upload node uniforms
	this->shader->setUniform("u_viewprojection", camera->viewprojection_matrix);
	this->shader->setUniform("u_camera_position", camera->eye);
	this->shader->setUniform("u_model", model);
	

	this->shader->setUniform("u_color", this->color);

	if (this->texture) {
		this->shader->setUniform("u_texture", this->texture);
	}
}

void StandardMaterial::render(Mesh* mesh, glm::mat4 model, Camera* camera)
{
	bool first_pass = true;
	if (mesh && this->shader)
	{
		// enable shader
		this->shader->enable();

		// Multi pass render
		int num_lights = Application::instance->light_list.size();
		for (int nlight = -1; nlight < num_lights; nlight++)
		{
			if (nlight == -1) { nlight++; } // hotfix

			// upload uniforms
			setUniforms(camera, model);

			// upload light uniforms
			if (!first_pass) {
				glBlendFunc(GL_SRC_ALPHA, GL_ONE);
				glDepthFunc(GL_LEQUAL);
			}
			this->shader->setUniform("u_ambient_light", Application::instance->ambient_light * (float)first_pass);

			if (num_lights > 0) {
				Light* light = Application::instance->light_list[nlight];
				light->setUniforms(this->shader, model);
			}
			else {
				// Set some uniforms in case there is no light
				this->shader->setUniform("u_light_intensity", 1.f);
				this->shader->setUniform("u_light_shininess", 1.f);
				this->shader->setUniform("u_light_color", glm::vec4(0.f));
			}

			// do the draw call
			mesh->render(GL_TRIANGLES);

			first_pass = false;
		}

		// disable shader
		this->shader->disable();
	}
}

void StandardMaterial::renderInMenu()
{
	if (ImGui::Checkbox("Show Normals", &this->show_normals)) {
		if (this->show_normals) {
			this->shader = this->normal_shader;
		}
		else {
			this->shader = this->base_shader;
		}
	}

	if (!this->show_normals) ImGui::ColorEdit3("Color", (float*)&this->color);
}

VolumeMaterial::~VolumeMaterial() { }

VolumeMaterial::VolumeMaterial(glm::vec4 color) {

	this->color = color;

	//Load shaders

	loadVDB("res/volumes/bunny_cloud.vdb"); 
	this->shader = Shader::Get("res/shaders/basic.vs", "res/shaders/volume_bunny.fs");

	//Deafult values

	this->density_mode = Bunny;

	this->absortion_coefitient = 0.9f; 
	step_length = 0.05f;

}

void VolumeMaterial::setUniforms(Camera* camera, glm::mat4 model)
{

	//upload node uniforms
	this->shader->setUniform("u_viewprojection", camera->viewprojection_matrix);
	this->shader->setUniform("u_camera_position", camera->eye);
	this->shader->setUniform("u_model", model);

	glm::mat4 inv_model = glm::inverse(model);
	glm::vec4 tmp = glm::vec4(camera->eye, 1.0f);
	tmp = inv_model * tmp;
	glm::vec3 camera_pos_local = (1.0f / tmp.w) * glm::vec3(tmp.x, tmp.y, tmp.z);

	this->shader->setUniform("u_local_camera_position", camera_pos_local);

	//printf("|%f, %f, %f|", this->color.x, this->color.y, this->color.z); 
	glm::vec4 color = glm::vec4(this->color.x, this->color.y, this->color.z, 1.0f); 
	this->shader->setUniform("u_color", color);

	// DEFINE BACKGROUND COLOR
	this->shader->setUniform("u_bg_color", Application::instance->background_color);


	glm::vec3 box_min = glm::vec3(-1, -1, -1);
	glm::vec3 box_max = glm::vec3(1, 1, 1);
	this->shader->setUniform("u_box_min", box_min);
	this->shader->setUniform("u_box_max", box_max);


	this->shader->setUniform("u_absortion_coef", absortion_coefitient);
	this->shader->setUniform("u_step_length", step_length);
	this->shader->setUniform("u_scale", scale);
	this->shader->setUniform("u_detail", detail);
	

	const unsigned int TEXTURE_SLOT = 0; 
	this->shader->setUniform("u_texture", this->texture, TEXTURE_SLOT);
	this->shader->setUniform("u_source_density", this->density_mode);
	

	
	//LIGHTS
	glm::vec4 color2 = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
	glm::vec3 light_pos_local = glm::vec3(1.5f, 1.5f, -1.5f);
	
	
	// FOR LOOP FOR LIGHT SOURCES
	//es queda sempre l'ultim valor
	float intensity = 1.111f;
	for (int l = 0; l < MAX_LIGHT; l++) {
		this->shader->setUniform("u_light_intensity", intensity);
		this->shader->setUniform("u_light_color", color2);
		this->shader->setUniform("u_light_pos_local", light_pos_local);
	}

}

void VolumeMaterial::render(Mesh* mesh, glm::mat4 model, Camera* camera)
{

	if (mesh == NULL || this->shader == NULL) {
		return; 
	}

	this->shader->enable();

	//glBlendFunc(GL_SRC_ALPHA, GL_ONE);
	glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
	/*
	source = new_color
	d = framebuffer

	source * 1 + framebuffer * (1-alfa)
		
	(1-alfa) = (1-( 1.0 - transmitansse)) = transmitansse

	...
	*/
	glEnable(GL_BLEND);

	glDepthFunc(GL_LEQUAL);

	setUniforms(camera, model);

	mesh->render(GL_TRIANGLES);


	glDisable(GL_BLEND);
	this->shader->disable();



}

void VolumeMaterial::renderInMenu()
{

	ImGui::ColorEdit3("Color"		, (float*)&this->color);
	ImGui::SliderFloat("Absortion coefitient"	, &this->absortion_coefitient, 0.0f, 4.0f);
	ImGui::SliderFloat("Step length", &this->step_length, 0.004f, 1.0f);
	ImGui::SliderFloat("Scale"		, &this->scale		, 0.001f, 4.5f);
	ImGui::SliderFloat("Detail"		, &this->detail		, 0.001f, 8.0f);

	ImGui::Combo("Density mode"		, (int*)&this->density_mode, "Homogeneus\0Noise\0Bunny\0"); 

}


void VolumeMaterial::loadVDB(std::string file_path)
{
	easyVDB::OpenVDBReader* vdbReader = new easyVDB::OpenVDBReader();
	vdbReader->read(file_path);

	// now, read the grid from the vdbReader and store the data in a 3D texture
	estimate3DTexture(vdbReader);
}

void VolumeMaterial::estimate3DTexture(easyVDB::OpenVDBReader* vdbReader)
{
	int resolution = 128;
	float radius = 2.0;

	int convertedGrids = 0;
	int convertedVoxels = 0;

	int totalGrids = vdbReader->gridsSize;
	int totalVoxels = totalGrids * pow(resolution, 3);

	float resolutionInv = 1.0f / resolution;
	int resolutionPow2 = pow(resolution, 2);
	int resolutionPow3 = pow(resolution, 3);

	// read all grids data and convert to texture
	for (unsigned int i = 0; i < totalGrids; i++) {
		easyVDB::Grid& grid = vdbReader->grids[i];
		float* data = new float[resolutionPow3];
		memset(data, 0, sizeof(float) * resolutionPow3);

		// Bbox
		easyVDB::Bbox bbox = easyVDB::Bbox();
		bbox = grid.getPreciseWorldBbox();
		glm::vec3 target = bbox.getCenter();
		glm::vec3 size = bbox.getSize();
		glm::vec3 step = size * resolutionInv;

		grid.transform->applyInverseTransformMap(step);
		target = target - (size * 0.5f);
		grid.transform->applyInverseTransformMap(target);
		target = target + (step * 0.5f);

		int x = 0;
		int y = 0;
		int z = 0;

		for (unsigned int j = 0; j < resolutionPow3; j++) {
			int baseX = x;
			int baseY = y;
			int baseZ = z;
			int baseIndex = baseX + baseY * resolution + baseZ * resolutionPow2;

			if (target.x >= 40 && target.y >= 40.33 && target.z >= 10.36) {
				int a = 0;
			}

			float value = grid.getValue(target);

			int cellBleed = radius;

			if (cellBleed) {
				for (int sx = -cellBleed; sx < cellBleed; sx++) {
					for (int sy = -cellBleed; sy < cellBleed; sy++) {
						for (int sz = -cellBleed; sz < cellBleed; sz++) {
							if (x + sx < 0.0 || x + sx >= resolution ||
								y + sy < 0.0 || y + sy >= resolution ||
								z + sz < 0.0 || z + sz >= resolution) {
								continue;
							}

							int targetIndex = baseIndex + sx + sy * resolution + sz * resolutionPow2;

							float offset = std::max(0.0, std::min(1.0, 1.0 - std::hypot(sx, sy, sz) / (radius / 2.0)));
							float dataValue = offset * value * 255.f;

							data[targetIndex] += dataValue;
							data[targetIndex] = std::min((float)data[targetIndex], 255.f);
						}
					}
				}
			}
			else {
				float dataValue = value * 255.f;

				data[baseIndex] += dataValue;
				data[baseIndex] = std::min((float)data[baseIndex], 255.f);
			}

			convertedVoxels++;

			if (z >= resolution) {
				break;
			}

			x++;
			target.x += step.x;

			if (x >= resolution) {
				x = 0;
				target.x -= step.x * resolution;

				y++;
				target.y += step.y;
			}

			if (y >= resolution) {
				y = 0;
				target.y -= step.y * resolution;

				z++;
				target.z += step.z;
			}

			// yield
		}

		// now we create the texture with the data
		// use this: https://www.khronos.org/opengl/wiki/OpenGL_Type
		// and this: https://registry.khronos.org/OpenGL-Refpages/gl4/html/glTexImage3D.xhtml
		this->texture = new Texture();
		this->texture->create3D(resolution, resolution, resolution, GL_RED, GL_FLOAT, false, data, GL_R8);
	}
}

