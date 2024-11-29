#include "application.h"

bool render_wireframe = false;
Camera* Application::camera = nullptr;

void Application::init(GLFWwindow* window)
{
    float to_rgb = 1.0f / 255.0f; 
    this->instance = this;
    glfwGetFramebufferSize(window, &this->window_width, &this->window_height);

    // OpenGL flags
    glEnable(GL_CULL_FACE); // render both sides of every triangle
    glEnable(GL_DEPTH_TEST); // check the occlusions using the Z buffer

    // Create camera
    this->camera = new Camera();
    this->camera->lookAt(glm::vec3(1.f, 1.5f, 4.f), glm::vec3(0.f, 0.0f, 0.f), glm::vec3(0.f, 1.f, 0.f));
    this->camera->setPerspective(60.f, this->window_width / (float)this->window_height, 0.1f, 500.f); // set the projection, we want to be perspective

    this->flag_grid = true;
    this->flag_wireframe = false;

    this->ambient_light = glm::vec4(0.75f, 0.75f, 0.75f, 1.f);
    bad_time = 0.0f; 

    /* ADD NODES TO THE SCENE */
    SceneNode* example = new SceneNode();
    example->mesh = Mesh::Get("res/meshes/cube.obj");
    bool use_normal_volume = false; 
    if (use_normal_volume) {
        example->material = new VolumeMaterial(glm::vec4(0, 13 * to_rgb, 0, 1));

        Light* sun = new Light(glm::vec3(1.5f, 1.5f, -1.5f), LIGHT_POINT, 12.8f, glm::vec4(1.0f, 143 * to_rgb, 65 * to_rgb, 1.0f));
        this->light_list.push_back(sun);
        this->node_list.push_back(sun);

        Light* moon = new Light(glm::vec3(1.5f, 2.0f, 0.0f), LIGHT_POINT, 3.141592f * 3.141592f, glm::vec4(50 * to_rgb, 100 * to_rgb, 200 * to_rgb, 1.0f));
        this->light_list.push_back(moon);
        this->node_list.push_back(moon);

    }
    else {
        glm::vec4 color = to_rgb * glm::vec4(240, 60, 14, 1);
        example->material = new IsosurfaceMaterial(color, 0.5);

        // simple white light
        Light* sun = new Light(glm::vec3(1.5f, 1.5f, -1.5f), LIGHT_POINT, 0.8, glm::vec4(1.0f, 1.0f, 1.0f, 1.0f));
        this->light_list.push_back(sun);
        this->node_list.push_back(sun);
    }

    this->node_list.push_back(example);

    //this->background_color = glm::vec3(219.0f / 255.0f, 237.0f / 255.0f, 242.0f / 255.0f);
    //this->background_color = glm::vec3(128.0f / 255.0f, 214.0f / 255.0f, 1.0f);
    this->background_color = glm::vec3(25 * to_rgb, 25 * to_rgb, 25 * to_rgb);




    this->speed = 0.1f;


}

void Application::update(float dt)
{
    // mouse update
    glm::vec2 delta = this->lastMousePosition - this->mousePosition;
    if (this->dragging) {
        this->camera->orbit(-delta.x * dt, delta.y * dt);
    }
    this->lastMousePosition = this->mousePosition;

    float time_mult = (2.0 * 355.0f / 113.0f); 
    //ImGui::SliderFloat("Light rotation speed", &this->speed, -1.0f, 1.0);

    float framerate = 60.0f; //fps
    bad_time += 1.0f / framerate;

    if (abs(this->speed) < 0.004f) {
        return; 
    }

    time_mult = time_mult * this->speed;

    float x_org = 1.5f;
    float z_org = -1.5f;
    float c = cos(bad_time * time_mult);
    float s = sin(bad_time * time_mult);

    float new_x = x_org * c - z_org * s;
    float new_z = z_org * c + x_org * s;

    light_list[0]->model[3][0] = new_x;
    light_list[0]->model[3][2] = new_z;



}

void Application::render()
{
    // set the clear color (the background color)
    //glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    
    glClearColor(background_color[0], background_color[1], background_color[2], 1.0f);

    // Clear the window and the depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // set flags
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    for (unsigned int i = 0; i < this->node_list.size(); i++)
    {
        this->node_list[i]->render(this->camera);

        if (this->flag_wireframe) this->node_list[i]->renderWireframe(this->camera);
    }

    // Draw the floor grid
    if (this->flag_grid) drawGrid();
}

void Application::renderGUI()
{
    if (ImGui::TreeNodeEx("Scene", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::ColorEdit3("Backgroiund color", (float*)&this->background_color);
        ImGui::ColorEdit3("Ambient light", (float*)&this->ambient_light);

        if (ImGui::TreeNode("Camera")) {
            this->camera->renderInMenu();
            ImGui::TreePop();
        }

        unsigned int count = 0;
        std::stringstream ss;
        for (auto& node : this->node_list) {
            ss << count;
            if (ImGui::TreeNode(node->name.c_str())) {
                node->renderInMenu();
                ImGui::TreePop();
            }
        }
        //ImGui::SliderFloat("Light rotation speed", &this->speed, -1.0f, 1.0);
        ImGui::DragFloat("Light rotation speed", &this->speed, 0.002f, -1.0f, 1.0);

        ImGui::TreePop();
    }
}

void Application::shutdown() { }

// keycodes: https://www.glfw.org/docs/3.3/group__keys.html
void Application::onKeyDown(int key, int scancode)
{
    switch (key) {
    case GLFW_KEY_ESCAPE: // quit
        close = true;
        break;
    case GLFW_KEY_R:
        Shader::ReloadAll();
        break;
    }
}

// keycodes: https://www.glfw.org/docs/3.3/group__keys.html
void Application::onKeyUp(int key, int scancode)
{
    switch (key) {
    case GLFW_KEY_T:
        std::cout << "T released" << std::endl;
        break;
    }
}

void Application::onRightMouseDown()
{
    this->dragging = true;
    this->lastMousePosition = this->mousePosition;
}

void Application::onRightMouseUp()
{
    this->dragging = false;
    this->lastMousePosition = this->mousePosition;
}

void Application::onLeftMouseDown()
{
    this->dragging = true;
    this->lastMousePosition = this->mousePosition;
}

void Application::onLeftMouseUp()
{
    this->dragging = false;
    this->lastMousePosition = this->mousePosition;
}

void Application::onMiddleMouseDown() { }

void Application::onMiddleMouseUp() { }

void Application::onMousePosition(double xpos, double ypos) { }

void Application::onScroll(double xOffset, double yOffset)
{
    int min = this->camera->min_fov;
    int max = this->camera->max_fov;

    if (yOffset < 0) {
        this->camera->fov += 4.f;
        if (this->camera->fov > max) {
            this->camera->fov = max;
        }
    }
    else {
        this->camera->fov -= 4.f;
        if (this->camera->fov < min) {
            this->camera->fov = min;
        }
    }
    this->camera->updateProjectionMatrix();
}