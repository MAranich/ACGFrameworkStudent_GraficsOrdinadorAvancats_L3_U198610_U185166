#include "material.h"

#include "application.h"

#include <istream>
#include <fstream>
#include <algorithm>


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

	// homogeneus
	shader_list.push_back(Shader::Get("res/shaders/basic.vs", "res/shaders/volume_homogeneus.fs"));
	//rand heterogeneus
	shader_list.push_back(Shader::Get("res/shaders/basic.vs", "res/shaders/volume_rand_hetero.fs"));
	//heterogeneus
	shader_list.push_back(Shader::Get("res/shaders/basic.vs", "res/shaders/volume_hetero.fs")); 

	this->current_shader = RandHeterogeneus; 
	this->shader = shader_list[(int)current_shader];

	this->absortion_coefitient = 0.3f; 
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
	glm::vec4 color_4 = glm::vec4(this->color.x, this->color.z, this->color.z, 1.0f); 
	this->shader->setUniform("u_color", color_4);

	glm::vec3 box_min = glm::vec3(-1, -1, -1);
	glm::vec3 box_max = glm::vec3(1, 1, 1);
	this->shader->setUniform("u_box_min", box_min);
	this->shader->setUniform("u_box_max", box_max);


	this->shader->setUniform("u_absortion_coef", absortion_coefitient);
	this->shader->setUniform("u_step_length", step_length);



}

void VolumeMaterial::render(Mesh* mesh, glm::mat4 model, Camera* camera)
{

	if (mesh == NULL) {
		return; 
	}
	if (this->shader == NULL) {
		return; 
	}

	this->shader = shader_list[(int)current_shader]; 

	this->shader->enable();

	//glBlendFunc(GL_SRC_ALPHA, GL_ONE);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);

	glDepthFunc(GL_LEQUAL);

	setUniforms(camera, model);

	mesh->render(GL_TRIANGLES);


	glDisable(GL_BLEND);
	this->shader->disable();

}

void VolumeMaterial::renderInMenu()
{

	ImGui::ColorEdit3("Color", (float*)&this->color);
	ImGui::SliderFloat("Absortion coefitient", &this->absortion_coefitient, 0.0f, 1.0f);
	ImGui::SliderFloat("Step length", &this->step_length, 0.001f, 1.0f);

	ImGui::Combo("Name", (int*)&this->current_shader, "Homogeneus\0Heterogeneus random\0Heterogeneus\0"); 

}
