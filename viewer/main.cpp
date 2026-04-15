/**
 * DistributionTool Viewer — OVITO-style 3D atom visualizer
 *
 * Usage:
 *   ./distool_viewer output.csv
 *   ./distool_viewer analyzed_dump.dump    (LAMMPS dump with dist_to_ref column)
 *
 * Controls:
 *   Left-drag   : orbit camera
 *   Middle-drag : pan
 *   Scroll      : zoom
 *   Click atom  : show info
 */

// ── OpenGL / window ──────────────────────────────────────────────────────────
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// ── Math ─────────────────────────────────────────────────────────────────────
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// ── UI ───────────────────────────────────────────────────────────────────────
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// ── Std ──────────────────────────────────────────────────────────────────────
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
#include <unistd.h>   // readlink, access

// ═══════════════════════════════════════════════════════════════════════════════
// Data structures
// ═══════════════════════════════════════════════════════════════════════════════

struct AtomRecord {
    int    id          = 0;
    int    type        = 1;
    float  x, y, z;
    float  dist_to_ref = 0.f;
    float  defect_prob = 0.f;
    int    defect_label = 0;  // 0=Lattice,1=Interstitial,2=VacAdj,3=TypeA,4=Unknown
};

static const char* defectName(int d) {
    switch (d) {
        case 0: return "Lattice";
        case 1: return "Interstitial";
        case 2: return "VacancyAdj";
        case 3: return "TypeA";
        default: return "Unknown";
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// File loading
// ═══════════════════════════════════════════════════════════════════════════════

static std::vector<AtomRecord> loadCSV(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open: " + path);

    std::string line;
    // Read header
    if (!std::getline(f, line))
        throw std::runtime_error("Empty file: " + path);

    // Parse column names from header (starts with #)
    std::string hdr = line;
    if (!hdr.empty() && hdr[0] == '#') hdr = hdr.substr(1);
    std::istringstream hs(hdr);
    std::vector<std::string> cols;
    { std::string t; while (hs >> t) cols.push_back(t); }

    auto colIdx = [&](std::initializer_list<const char*> names) -> int {
        for (auto* n : names)
            for (int i = 0; i < (int)cols.size(); ++i)
                if (cols[i] == n) return i;
        return -1;
    };

    int id_c  = colIdx({"id"});
    int ty_c  = colIdx({"type"});
    int x_c   = colIdx({"x"});
    int y_c   = colIdx({"y"});
    int z_c   = colIdx({"z"});
    int dr_c  = colIdx({"dist_to_ref"});
    int dp_c  = colIdx({"defect_prob"});
    int dl_c  = colIdx({"defect_label","defect_type"});

    if (x_c < 0 || y_c < 0 || z_c < 0 || dr_c < 0)
        throw std::runtime_error("CSV missing required columns (x y z dist_to_ref)");

    std::vector<AtomRecord> atoms;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::vector<std::string> v;
        { std::string t; while (ss >> t) v.push_back(t); }
        if (v.empty()) continue;

        if ((int)v.size() <= std::max({x_c, y_c, z_c})) continue; // malformed row

        AtomRecord a;
        auto get = [&](int c) -> const std::string& {
            static const std::string empty;
            return (c >= 0 && c < (int)v.size()) ? v[c] : empty;
        };
        try {
            if (id_c >= 0) a.id   = std::stoi(get(id_c));
            if (ty_c >= 0) a.type = std::stoi(get(ty_c));
            a.x = std::stof(v[x_c]);
            a.y = std::stof(v[y_c]);
            a.z = std::stof(v[z_c]);
            if (dr_c >= 0 && dr_c < (int)v.size()) a.dist_to_ref = std::stof(v[dr_c]);
            if (dp_c >= 0 && dp_c < (int)v.size()) a.defect_prob = std::stof(v[dp_c]);
            if (dl_c >= 0 && dl_c < (int)v.size()) {
                const auto& ds = v[dl_c];
                if (!ds.empty() && (std::isdigit((unsigned char)ds[0]) || ds[0] == '-'))
                    a.defect_label = std::clamp(std::stoi(ds), 0, 4);
                else {
                    if      (ds == "Lattice")      a.defect_label = 0;
                    else if (ds == "Interstitial") a.defect_label = 1;
                    else if (ds == "VacancyAdj")   a.defect_label = 2;
                    else if (ds == "TypeA")        a.defect_label = 3;
                    else                           a.defect_label = 4;
                }
            }
        } catch (...) { continue; }  // skip malformed rows silently
        atoms.push_back(a);
    }
    return atoms;
}

// Load analyzed LAMMPS dump (ITEM: ATOMS id type x y z dist_to_ref defect_prob defect_label)
static std::vector<AtomRecord> loadDump(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open: " + path);

    std::string line;
    std::vector<AtomRecord> atoms;

    while (std::getline(f, line)) {
        if (line.find("ITEM: ATOMS") == std::string::npos) continue;

        // Parse column names
        std::istringstream hs(line);
        std::vector<std::string> cols;
        { std::string t; while (hs >> t) cols.push_back(t); }
        // Remove "ITEM:" and "ATOMS"
        if (cols.size() >= 2) cols.erase(cols.begin(), cols.begin()+2);

        auto colIdx = [&](std::initializer_list<const char*> names) -> int {
            for (auto* n : names)
                for (int i = 0; i < (int)cols.size(); ++i)
                    if (cols[i] == n) return i;
            return -1;
        };
        int id_c = colIdx({"id"});
        int ty_c = colIdx({"type"});
        int x_c  = colIdx({"x","xu","xs"});
        int y_c  = colIdx({"y","yu","ys"});
        int z_c  = colIdx({"z","zu","zs"});
        int dr_c = colIdx({"dist_to_ref"});
        int dp_c = colIdx({"defect_prob"});
        int dl_c = colIdx({"defect_label"});

        if (x_c < 0 || y_c < 0 || z_c < 0)
            throw std::runtime_error("Dump missing x/y/z columns");

        while (std::getline(f, line)) {
            if (line.find("ITEM:") != std::string::npos) break;
            if (line.empty()) continue;
            std::istringstream ss(line);
            std::vector<std::string> v;
            { std::string t; while (ss >> t) v.push_back(t); }
            if (v.empty()) continue;

            AtomRecord a;
            auto get = [&](int c) -> const std::string& {
                static const std::string empty;
                return (c >= 0 && c < (int)v.size()) ? v[c] : empty;
            };
            if (id_c >= 0) a.id   = std::stoi(get(id_c));
            if (ty_c >= 0) a.type = std::stoi(get(ty_c));
            a.x = std::stof(get(x_c));
            a.y = std::stof(get(y_c));
            a.z = std::stof(get(z_c));
            if (dr_c >= 0) a.dist_to_ref = std::stof(get(dr_c));
            if (dp_c >= 0) a.defect_prob = std::stof(get(dp_c));
            if (dl_c >= 0) a.defect_label = std::stoi(get(dl_c));
            atoms.push_back(a);
        }
        break; // Only first frame
    }
    return atoms;
}

// ═══════════════════════════════════════════════════════════════════════════════
// Camera (arcball orbit)
// ═══════════════════════════════════════════════════════════════════════════════

struct Camera {
    glm::vec3 target   = {0, 0, 0};
    float     distance = 100.f;
    float     theta    = 0.f;                    // azimuth (radians)
    float     phi      = glm::radians(25.f);     // elevation (radians)
    float     fov      = glm::radians(45.f);

    glm::vec3 position() const {
        return target + distance * glm::vec3(
            std::cos(phi) * std::sin(theta),
            std::sin(phi),
            std::cos(phi) * std::cos(theta));
    }

    glm::mat4 view() const {
        return glm::lookAt(position(), target, glm::vec3(0, 1, 0));
    }

    glm::mat4 proj(float aspect) const {
        return glm::perspective(fov, aspect, distance * 0.001f, distance * 10.f);
    }

    void orbit(float dtheta, float dphi) {
        theta += dtheta;
        phi = glm::clamp(phi + dphi,
                         glm::radians(-89.f), glm::radians(89.f));
    }

    void pan(float dx, float dy) {
        glm::vec3 right = glm::normalize(
            glm::cross(target - position(), glm::vec3(0, 1, 0)));
        glm::vec3 up = glm::vec3(0, 1, 0);
        target += right * dx + up * dy;
    }

    void zoom(float factor) {
        distance = std::max(0.1f, distance * factor);
    }
};

// ═══════════════════════════════════════════════════════════════════════════════
// Shader helpers
// ═══════════════════════════════════════════════════════════════════════════════

static GLuint compileShader(GLenum type, const char* src) {
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char buf[1024]; glGetShaderInfoLog(s, 1024, nullptr, buf);
        throw std::runtime_error(std::string("Shader compile error: ") + buf);
    }
    return s;
}

static GLuint linkProgram(const char* vert_src, const char* frag_src) {
    GLuint v = compileShader(GL_VERTEX_SHADER,   vert_src);
    GLuint f = compileShader(GL_FRAGMENT_SHADER, frag_src);
    GLuint p = glCreateProgram();
    glAttachShader(p, v); glAttachShader(p, f);
    glLinkProgram(p);
    glDeleteShader(v); glDeleteShader(f);
    GLint ok; glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) {
        char buf[1024]; glGetProgramInfoLog(p, 1024, nullptr, buf);
        throw std::runtime_error(std::string("Program link error: ") + buf);
    }
    return p;
}

// ── Sphere impostor shaders ───────────────────────────────────────────────────

static const char* VERT_SRC = R"glsl(
#version 330 core

// Per-vertex: unit quad corner
layout(location = 0) in vec2 a_corner;

// Per-instance
layout(location = 1) in vec3  a_center;
layout(location = 2) in float a_radius;
layout(location = 3) in vec4  a_color;

out vec3  v_center_view;
out float v_radius;
out vec4  v_color;
out vec2  v_corner;

uniform mat4 u_mv;
uniform mat4 u_proj;

void main() {
    v_center_view = vec3(u_mv * vec4(a_center, 1.0));
    v_radius      = a_radius;
    v_color       = a_color;
    v_corner      = a_corner;

    // Billboard: expand quad by radius in view space
    vec3 pos_view = v_center_view + vec3(a_corner * a_radius, 0.0);
    gl_Position   = u_proj * vec4(pos_view, 1.0);
}
)glsl";

static const char* FRAG_SRC = R"glsl(
#version 330 core

in vec3  v_center_view;
in float v_radius;
in vec4  v_color;
in vec2  v_corner;

out vec4 frag_color;

uniform mat4 u_proj;

void main() {
    float r2 = dot(v_corner, v_corner);
    if (r2 > 1.0) discard;

    // Surface normal in view space
    vec3 normal = vec3(v_corner, sqrt(1.0 - r2));

    // Hit point in view space
    vec3 hit = v_center_view + normal * v_radius;

    // Phong lighting (fixed light in view space)
    vec3  L        = normalize(vec3(0.6, 1.0, 0.8));
    float diffuse  = max(dot(normal, L), 0.0);
    float specular = pow(max(dot(reflect(-L, normal), vec3(0,0,1)), 0.0), 32.0);
    float light    = 0.25 + 0.65 * diffuse + 0.10 * specular;

    frag_color = vec4(v_color.rgb * light, v_color.a);

    // Correct depth so spheres occlude each other properly
    vec4 clip    = u_proj * vec4(hit, 1.0);
    gl_FragDepth = (clip.z / clip.w + 1.0) * 0.5;
}
)glsl";

// ── Picking shader (flat color per instance for mouse picking) ────────────────

static const char* PICK_VERT_SRC = R"glsl(
#version 330 core
layout(location = 0) in vec2  a_corner;
layout(location = 1) in vec3  a_center;
layout(location = 2) in float a_radius;
layout(location = 3) in vec4  a_pick_color; // encodes atom index as RGBA uint8
out vec3  v_center_view;
out float v_radius;
out vec4  v_pick_color;
out vec2  v_corner;
uniform mat4 u_mv;
uniform mat4 u_proj;
void main() {
    v_center_view = vec3(u_mv * vec4(a_center, 1.0));
    v_radius      = a_radius;
    v_pick_color  = a_pick_color;
    v_corner      = a_corner;
    vec3 pos_view = v_center_view + vec3(a_corner * a_radius, 0.0);
    gl_Position   = u_proj * vec4(pos_view, 1.0);
}
)glsl";

static const char* PICK_FRAG_SRC = R"glsl(
#version 330 core
in vec3  v_center_view;
in float v_radius;
in vec4  v_pick_color;
in vec2  v_corner;
out vec4 frag_color;
uniform mat4 u_proj;
void main() {
    float r2 = dot(v_corner, v_corner);
    if (r2 > 1.0) discard;
    vec3 normal = vec3(v_corner, sqrt(1.0 - r2));
    vec3 hit    = v_center_view + normal * v_radius;
    vec4 clip   = u_proj * vec4(hit, 1.0);
    gl_FragDepth = (clip.z / clip.w + 1.0) * 0.5;
    frag_color  = v_pick_color;
}
)glsl";

// ═══════════════════════════════════════════════════════════════════════════════
// GPU instance data
// ═══════════════════════════════════════════════════════════════════════════════

struct GpuAtom {
    float x, y, z;       // center
    float radius;
    float r, g, b, a;    // color
};

// ═══════════════════════════════════════════════════════════════════════════════
// Color mapping
// ═══════════════════════════════════════════════════════════════════════════════

enum class ColorMode { DistToRef = 0, DefectType, AtomType, DefectProb };

static glm::vec4 plasma(float t) {
    // Plasma colormap approximation
    t = glm::clamp(t, 0.f, 1.f);
    float r = 0.050f + t * (1.700f - t * 0.700f);
    float g = 0.030f + t * (0.500f - t * 1.100f) + t*t * 0.700f;
    float b = 0.530f + t * (1.000f - t * 2.300f) + t*t * 1.400f;
    return {glm::clamp(r,0.f,1.f), glm::clamp(g,0.f,1.f), glm::clamp(b,0.f,1.f), 1.f};
}

static const glm::vec4 DEFECT_COLORS[5] = {
    {0.45f, 0.75f, 0.45f, 1.f},  // Lattice — green
    {0.95f, 0.30f, 0.25f, 1.f},  // Interstitial — red
    {0.25f, 0.55f, 0.95f, 1.f},  // VacancyAdj — blue
    {0.95f, 0.75f, 0.15f, 1.f},  // TypeA — yellow
    {0.70f, 0.35f, 0.90f, 1.f},  // Unknown — purple
};

static const glm::vec4 TYPE_COLORS[] = {
    {0.90f, 0.40f, 0.10f, 1.f},
    {0.10f, 0.55f, 0.90f, 1.f},
    {0.20f, 0.80f, 0.40f, 1.f},
    {0.95f, 0.80f, 0.10f, 1.f},
    {0.80f, 0.20f, 0.80f, 1.f},
    {0.10f, 0.80f, 0.80f, 1.f},
    {0.80f, 0.50f, 0.10f, 1.f},
    {0.50f, 0.50f, 0.50f, 1.f},
};

// ═══════════════════════════════════════════════════════════════════════════════
// Renderer
// ═══════════════════════════════════════════════════════════════════════════════

class Renderer {
public:
    GLuint prog      = 0;
    GLuint pick_prog = 0;
    GLuint vao = 0, quad_vbo = 0, inst_vbo = 0;
    GLuint pick_vao = 0, pick_inst_vbo = 0;
    GLuint pick_fbo = 0, pick_tex = 0, pick_rbo = 0;
    int    inst_count = 0;

    void init() {
        prog      = linkProgram(VERT_SRC,      FRAG_SRC);
        pick_prog = linkProgram(PICK_VERT_SRC, PICK_FRAG_SRC);

        // Unit quad
        float quad[] = { -1,-1,  1,-1,  -1,1,  1,1 };
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        glGenBuffers(1, &quad_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, quad_vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quad), quad, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, nullptr);

        glGenBuffers(1, &inst_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, inst_vbo);
        // location 1: center (3 floats), location 2: radius (1), location 3: color (4)
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GpuAtom), (void*)0);
        glVertexAttribDivisor(1, 1);
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(GpuAtom), (void*)(3*4));
        glVertexAttribDivisor(2, 1);
        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(GpuAtom), (void*)(4*4));
        glVertexAttribDivisor(3, 1);

        // Pick VAO (same quad, different instance buffer)
        glGenVertexArrays(1, &pick_vao);
        glBindVertexArray(pick_vao);
        glBindBuffer(GL_ARRAY_BUFFER, quad_vbo);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, nullptr);

        glGenBuffers(1, &pick_inst_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, pick_inst_vbo);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GpuAtom), (void*)0);
        glVertexAttribDivisor(1, 1);
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(GpuAtom), (void*)(3*4));
        glVertexAttribDivisor(2, 1);
        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(GpuAtom), (void*)(4*4));
        glVertexAttribDivisor(3, 1);

        glBindVertexArray(0);
    }

    void initPickFBO(int w, int h) {
        if (pick_fbo) {
            glDeleteFramebuffers(1, &pick_fbo);
            glDeleteTextures(1, &pick_tex);
            glDeleteRenderbuffers(1, &pick_rbo);
            pick_fbo = pick_tex = pick_rbo = 0;
        }
        glGenFramebuffers(1, &pick_fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, pick_fbo);

        glGenTextures(1, &pick_tex);
        glBindTexture(GL_TEXTURE_2D, pick_tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, pick_tex, 0);

        glGenRenderbuffers(1, &pick_rbo);
        glBindRenderbuffer(GL_RENDERBUFFER, pick_rbo);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, pick_rbo);

        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
            std::cerr << "Warning: pick FBO incomplete\n";

        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    void upload(const std::vector<GpuAtom>& data) {
        inst_count = (int)data.size();
        glBindBuffer(GL_ARRAY_BUFFER, inst_vbo);
        glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(GpuAtom), data.data(), GL_DYNAMIC_DRAW);
    }

    void uploadPick(const std::vector<GpuAtom>& data) {
        glBindBuffer(GL_ARRAY_BUFFER, pick_inst_vbo);
        glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(GpuAtom), data.data(), GL_DYNAMIC_DRAW);
    }

    void draw(const glm::mat4& mv, const glm::mat4& proj) {
        glUseProgram(prog);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_mv"),   1, GL_FALSE, glm::value_ptr(mv));
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_proj"), 1, GL_FALSE, glm::value_ptr(proj));
        glBindVertexArray(vao);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, inst_count);
        glBindVertexArray(0);
    }

    // Returns atom index (0-based in filtered list) or -1 if background
    int pick(int mouse_x, int mouse_y, int fb_h,
             const glm::mat4& mv, const glm::mat4& proj) {
        glBindFramebuffer(GL_FRAMEBUFFER, pick_fbo);
        glClearColor(1,1,1,1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(pick_prog);
        glUniformMatrix4fv(glGetUniformLocation(pick_prog, "u_mv"),   1, GL_FALSE, glm::value_ptr(mv));
        glUniformMatrix4fv(glGetUniformLocation(pick_prog, "u_proj"), 1, GL_FALSE, glm::value_ptr(proj));
        glBindVertexArray(pick_vao);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, inst_count);
        glBindVertexArray(0);

        unsigned char px[4];
        int y_flip = fb_h - 1 - mouse_y;
        glReadPixels(mouse_x, y_flip, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, px);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        if (px[0] == 255 && px[1] == 255 && px[2] == 255) return -1;
        int idx = px[0] | (px[1] << 8) | (px[2] << 16);
        return idx;
    }
};

// ═══════════════════════════════════════════════════════════════════════════════
// Analysis pipeline state
// ═══════════════════════════════════════════════════════════════════════════════

struct AnalysisState {
    // Paths
    char ref_path[1024] = {};
    char dmg_path[1024] = {};
    char out_path[512]  = "output.csv";

    // SOAP
    int   n_max = 9;
    int   l_max = 9;
    float r_cut = 5.0f;

    // Classification
    float threshold  = 0.15f;
    float vac_dist   = -1.0f;   // -1 = auto

    // Output options
    bool do_pca      = false;
    bool do_hist     = false;
    bool auto_load   = true;    // switch to viewer tab when done

    // Runtime state (accessed from worker thread)
    std::atomic<bool>   running{false};
    std::atomic<bool>   done{false};
    bool                success       = false;
    bool                result_loaded = false;
    bool                scroll_log    = false;
    std::string         result_csv;
    std::vector<std::string> log;
    std::mutex               log_mutex;
    std::thread              worker;
};

// ═══════════════════════════════════════════════════════════════════════════════
// Application state
// ═══════════════════════════════════════════════════════════════════════════════

struct App {
    // Data
    std::vector<AtomRecord> all_atoms;
    std::vector<int>        visible_indices;   // indices into all_atoms

    // Filter state
    float  dist_min = 0.f, dist_max = 1.f;
    float  filter_lo = 0.f, filter_hi = 1.f;
    bool   defect_filter[5]  = {true,true,true,true,true};
    std::vector<int>  atom_types;
    std::unordered_map<int,bool> type_visible;
    ColorMode color_mode = ColorMode::DistToRef;
    float  atom_radius   = 1.0f;
    float  opacity       = 1.0f;

    // Camera
    Camera camera;

    // Mouse state
    bool   mouse_left   = false;
    bool   mouse_mid    = false;
    double last_mx = 0, last_my = 0;
    int    win_w = 1400, win_h = 900;

    // Selection
    int    selected_idx = -1;   // index in visible_indices

    // Stats
    std::string filename;

    // Sub-systems
    AnalysisState analysis;
    int active_tab = 0;   // 0=Analizar, 1=Visualizar
};

// ═══════════════════════════════════════════════════════════════════════════════
// Build GPU data from filter state
// ═══════════════════════════════════════════════════════════════════════════════

static void buildGpuData(App& app, Renderer& rend) {
    app.visible_indices.clear();

    // Gather filtered atoms
    for (int i = 0; i < (int)app.all_atoms.size(); ++i) {
        const auto& a = app.all_atoms[i];
        if (a.dist_to_ref < app.filter_lo || a.dist_to_ref > app.filter_hi) continue;
        if (!app.defect_filter[std::min(a.defect_label, 4)])              continue;
        auto it = app.type_visible.find(a.type);
        if (it != app.type_visible.end() && !it->second)                  continue;
        app.visible_indices.push_back(i);
    }

    float drange = app.dist_max - app.dist_min;
    if (drange < 1e-9f) drange = 1.f;

    std::vector<GpuAtom> data, pick_data;
    data.reserve(app.visible_indices.size());
    pick_data.reserve(app.visible_indices.size());

    for (int vi = 0; vi < (int)app.visible_indices.size(); ++vi) {
        const auto& a = app.all_atoms[app.visible_indices[vi]];

        glm::vec4 col;
        switch (app.color_mode) {
            case ColorMode::DistToRef: {
                float t = (a.dist_to_ref - app.dist_min) / drange;
                col = plasma(t);
                break;
            }
            case ColorMode::DefectType:
                col = DEFECT_COLORS[std::min(a.defect_label, 4)];
                break;
            case ColorMode::AtomType: {
                int ti = (a.type - 1) % 8;
                col = TYPE_COLORS[ti < 0 ? 0 : ti];
                break;
            }
            case ColorMode::DefectProb: {
                col = plasma(a.defect_prob);
                break;
            }
        }
        col.a = app.opacity;

        data.push_back({a.x, a.y, a.z, app.atom_radius,
                        col.r, col.g, col.b, col.a});

        // Encode index as RGB (supports up to 16M atoms)
        unsigned char pr = vi & 0xFF, pg = (vi>>8)&0xFF, pb = (vi>>16)&0xFF;
        pick_data.push_back({a.x, a.y, a.z, app.atom_radius,
                             pr/255.f, pg/255.f, pb/255.f, 1.f});
    }

    rend.upload(data);
    rend.uploadPick(pick_data);
    app.selected_idx = -1;
}

// ═══════════════════════════════════════════════════════════════════════════════
// Center camera on loaded data
// ═══════════════════════════════════════════════════════════════════════════════

static void centerCamera(App& app) {
    if (app.all_atoms.empty()) return;
    glm::vec3 mn(1e30f), mx(-1e30f);
    for (const auto& a : app.all_atoms) {
        mn = glm::min(mn, glm::vec3(a.x, a.y, a.z));
        mx = glm::max(mx, glm::vec3(a.x, a.y, a.z));
    }
    app.camera.target   = (mn + mx) * 0.5f;
    app.camera.distance = glm::length(mx - mn) * 0.7f;
}

// ═══════════════════════════════════════════════════════════════════════════════
// Analysis helpers
// ═══════════════════════════════════════════════════════════════════════════════

// Opens a native file-picker via zenity (Linux). Falls back to empty string
// if zenity is not available — user can type the path manually.
static std::string browseFile(const char* title) {
    std::string cmd = std::string("zenity --file-selection --title=\"") + title + "\" 2>/dev/null";
    FILE* f = popen(cmd.c_str(), "r");
    if (!f) return "";
    char buf[4096] = {};
    if (fgets(buf, sizeof(buf), f)) { /* got path */ }
    pclose(f);
    std::string s(buf);
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    return s;
}

// Locate the distool executable relative to this viewer or in PATH.
static std::string findDistool() {
    char exe[4096] = {};
    ssize_t len = readlink("/proc/self/exe", exe, sizeof(exe) - 1);
    if (len > 0) {
        std::string dir(exe, (size_t)len);
        dir = dir.substr(0, dir.rfind('/'));
        for (const char* rel : {
                "/distool",
                "/../build_fresh/distool",
                "/../build/distool",
                "/../build_audit/distool",
                "/../build_viewer/../build_fresh/distool"}) {
            std::string p = dir + rel;
            if (access(p.c_str(), X_OK) == 0) return p;
        }
    }
    // Try PATH
    FILE* f = popen("which distool 2>/dev/null", "r");
    if (f) {
        char buf[512] = {};
        if (fgets(buf, sizeof(buf), f)) { /* ok */ }
        pclose(f);
        std::string s(buf);
        while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
        if (!s.empty()) return s;
    }
    return "";
}

// Runs in a background thread — builds and executes the distool command,
// piping stdout+stderr to AnalysisState::log line by line.
static void runAnalysisThread(AnalysisState* st) {
    std::string dtool = findDistool();
    if (dtool.empty()) {
        std::lock_guard<std::mutex> lk(st->log_mutex);
        st->log.push_back("[ERROR] distool no encontrado.");
        st->log.push_back("  Compilalo con: cmake --build build_fresh -j$(nproc)");
        st->success = false; st->done = true; st->running = false;
        return;
    }

    // Build command string
    auto fstr = [](float v) {
        char buf[32]; snprintf(buf, sizeof(buf), "%.4f", (double)v); return std::string(buf);
    };
    std::string cmd;
    cmd += "\"" + dtool + "\"";
    cmd += " \"" + std::string(st->ref_path) + "\"";
    cmd += " \"" + std::string(st->dmg_path) + "\"";
    cmd += " --n-max "    + std::to_string(st->n_max);
    cmd += " --l-max "    + std::to_string(st->l_max);
    cmd += " --r-cut "    + fstr(st->r_cut);
    cmd += " --threshold "+ fstr(st->threshold);
    if (st->vac_dist > 0)
        cmd += " --vac-dist " + fstr(st->vac_dist);
    cmd += " --output \"" + std::string(st->out_path) + "\"";
    if (st->do_pca)  cmd += " --pca";
    if (st->do_hist) cmd += " --hist";
    cmd += " 2>&1";

    {
        std::lock_guard<std::mutex> lk(st->log_mutex);
        st->log.push_back("$ " + cmd);
        st->log.push_back("");
        st->scroll_log = true;
    }

    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        std::lock_guard<std::mutex> lk(st->log_mutex);
        st->log.push_back("[ERROR] popen falló.");
        st->success = false; st->done = true; st->running = false;
        return;
    }

    char buf[512];
    while (fgets(buf, sizeof(buf), pipe)) {
        std::string line(buf);
        while (!line.empty() && (line.back()=='\n'||line.back()=='\r')) line.pop_back();
        std::lock_guard<std::mutex> lk(st->log_mutex);
        st->log.push_back(std::move(line));
        st->scroll_log = true;
    }

    int ret = pclose(pipe);
    st->success    = (ret == 0);
    st->result_csv = std::string(st->out_path);
    st->done       = true;
    st->running    = false;
}

// ═══════════════════════════════════════════════════════════════════════════════
// GLFW callbacks
// ═══════════════════════════════════════════════════════════════════════════════

static App*      g_app  = nullptr;
static Renderer* g_rend = nullptr;

static void cbScroll(GLFWwindow*, double, double dy) {
    if (!g_app || ImGui::GetIO().WantCaptureMouse) return;
    g_app->camera.zoom(dy > 0 ? 0.88f : 1.14f);
}

static void cbMouseBtn(GLFWwindow* w, int btn, int action, int) {
    if (!g_app) return;
    double mx, my;
    glfwGetCursorPos(w, &mx, &my);
    g_app->last_mx = mx;
    g_app->last_my = my;

    if (ImGui::GetIO().WantCaptureMouse) return;

    if (btn == GLFW_MOUSE_BUTTON_LEFT) {
        g_app->mouse_left = (action == GLFW_PRESS);
        if (action == GLFW_PRESS) {
            // Pick — convert window coords → framebuffer coords (HiDPI safe)
            int fbw, fbh, winw, winh;
            glfwGetFramebufferSize(w, &fbw, &fbh);
            glfwGetWindowSize(w, &winw, &winh);
            int fx = (int)(mx * (double)fbw / winw);
            int fy = (int)(my * (double)fbh / winh);
            auto mv   = g_app->camera.view();
            auto proj = g_app->camera.proj((float)fbw / (float)fbh);
            int vi = g_rend->pick(fx, fy, fbh, mv, proj);
            g_app->selected_idx = (vi >= 0 && vi < (int)g_app->visible_indices.size()) ? vi : -1;
        }
    }
    if (btn == GLFW_MOUSE_BUTTON_MIDDLE)
        g_app->mouse_mid = (action == GLFW_PRESS);
}

static void cbCursorPos(GLFWwindow*, double mx, double my) {
    if (!g_app) return;
    float dx = (float)(mx - g_app->last_mx);
    float dy = (float)(my - g_app->last_my);
    // Always update — prevents jump when cursor re-enters the 3D viewport
    g_app->last_mx = mx;
    g_app->last_my = my;

    if (ImGui::GetIO().WantCaptureMouse) return;

    if (g_app->mouse_left)
        g_app->camera.orbit(dx * 0.005f, -dy * 0.005f);
    else if (g_app->mouse_mid)
        g_app->camera.pan(-dx * g_app->camera.distance * 0.001f,
                           dy * g_app->camera.distance * 0.001f);
}

static void cbFramebuffer(GLFWwindow*, int w, int h) {
    if (g_rend) g_rend->initPickFBO(w, h);
    glViewport(0, 0, w, h);
}

// ═══════════════════════════════════════════════════════════════════════════════
// UI panel
// ═══════════════════════════════════════════════════════════════════════════════

// Load a CSV (or analyzed dump) into the viewer, updating all state.
static void loadIntoViewer(App& app, Renderer& rend, const std::string& path) {
    bool is_dump = path.find(".dump") != std::string::npos
                || path.rfind("dump", 0) != std::string::npos;
    std::vector<AtomRecord> atoms = is_dump ? loadDump(path) : loadCSV(path);
    if (atoms.empty()) return;

    app.all_atoms = std::move(atoms);
    app.filename  = path;

    app.dist_min = app.dist_max = app.all_atoms[0].dist_to_ref;
    for (const auto& a : app.all_atoms) {
        app.dist_min = std::min(app.dist_min, a.dist_to_ref);
        app.dist_max = std::max(app.dist_max, a.dist_to_ref);
    }
    app.filter_lo = app.dist_min;
    app.filter_hi = app.dist_max;

    app.atom_types.clear();
    app.type_visible.clear();
    for (const auto& a : app.all_atoms)
        if (app.type_visible.find(a.type) == app.type_visible.end()) {
            app.atom_types.push_back(a.type);
            app.type_visible[a.type] = true;
        }
    std::sort(app.atom_types.begin(), app.atom_types.end());

    for (int i = 0; i < 5; ++i) app.defect_filter[i] = true;

    centerCamera(app);
    buildGpuData(app, rend);
}


// ── Analysis tab ──────────────────────────────────────────────────────────────
static void drawAnalysisUI(App& app, Renderer& rend) {
    AnalysisState& st = app.analysis;
    const float W = -52;  // input width leaving room for "..." button

    // ── Reference dump ────────────────────────────────────────────────
    ImGui::TextDisabled("REFERENCIA (.dump cristal perfecto)");
    ImGui::SetNextItemWidth(W);
    ImGui::InputText("##ref", st.ref_path, sizeof(st.ref_path));
    ImGui::SameLine();
    if (ImGui::Button("...##r")) {
        auto p = browseFile("Seleccionar referencia.dump");
        if (!p.empty()) { strncpy(st.ref_path, p.c_str(), sizeof(st.ref_path)-1); }
    }

    // ── Damaged dump ──────────────────────────────────────────────────
    ImGui::TextDisabled("DAÑADO (.dump a analizar)");
    ImGui::SetNextItemWidth(W);
    ImGui::InputText("##dmg", st.dmg_path, sizeof(st.dmg_path));
    ImGui::SameLine();
    if (ImGui::Button("...##d")) {
        auto p = browseFile("Seleccionar damaged.dump");
        if (!p.empty()) { strncpy(st.dmg_path, p.c_str(), sizeof(st.dmg_path)-1); }
    }

    // ── Output ────────────────────────────────────────────────────────
    ImGui::TextDisabled("ARCHIVO DE SALIDA");
    ImGui::SetNextItemWidth(-1);
    ImGui::InputText("##out", st.out_path, sizeof(st.out_path));

    ImGui::Separator();

    // ── SOAP ─────────────────────────────────────────────────────────
    if (ImGui::CollapsingHeader("Configuración SOAP")) {
        ImGui::SliderInt("n_max",     &st.n_max, 1, 15);
        ImGui::SliderInt("l_max",     &st.l_max, 1, 15);
        ImGui::SliderFloat("r_cut (Å)", &st.r_cut, 1.f, 12.f, "%.1f");
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Radio de corte para descriptores SOAP");
        ImGui::TextDisabled("DV size: %d", (st.n_max*(st.n_max+1)/2) * (st.l_max+1));
    }

    // ── Defect detection ──────────────────────────────────────────────
    if (ImGui::CollapsingHeader("Detección de defectos")) {
        ImGui::SliderFloat("threshold", &st.threshold, 0.0f, 1.0f, "%.3f");
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Distancia mínima al DV de referencia para clasificar como defecto");
        ImGui::SliderFloat("vac_dist (Å)", &st.vac_dist, -1.f, 6.f, st.vac_dist < 0 ? "auto" : "%.2f");
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Threshold de detección de vacancias (-1 = auto: r_cut × 0.4)");
        ImGui::Checkbox("PCA",        &st.do_pca);
        ImGui::SameLine();
        ImGui::Checkbox("Histograma", &st.do_hist);
    }

    ImGui::Separator();
    ImGui::Checkbox("Cargar resultado al terminar", &st.auto_load);
    ImGui::Spacing();

    // ── RUN button ────────────────────────────────────────────────────
    bool can_run = !st.running
                   && st.ref_path[0] != '\0'
                   && st.dmg_path[0] != '\0'
                   && st.out_path[0] != '\0';

    if (!can_run) ImGui::BeginDisabled();
    ImGui::PushStyleColor(ImGuiCol_Button,        {0.12f, 0.50f, 0.12f, 1.f});
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, {0.18f, 0.68f, 0.18f, 1.f});
    ImGui::PushStyleColor(ImGuiCol_ButtonActive,  {0.08f, 0.35f, 0.08f, 1.f});
    bool clicked = ImGui::Button("  EJECUTAR ANÁLISIS  ", {-1, 36});
    ImGui::PopStyleColor(3);
    if (!can_run) ImGui::EndDisabled();

    if (clicked && can_run) {
        {
            std::lock_guard<std::mutex> lk(st.log_mutex);
            st.log.clear();
        }
        st.done = false; st.success = false; st.result_loaded = false;
        st.running = true;
        if (st.worker.joinable()) st.worker.join();
        st.worker = std::thread(runAnalysisThread, &st);
    }

    // ── Status ───────────────────────────────────────────────────────
    ImGui::Spacing();
    if (st.running) {
        ImGui::TextColored({1.f, 0.8f, 0.f, 1.f}, "  Procesando...");
    } else if (st.done) {
        if (st.success) {
            ImGui::TextColored({0.3f, 1.f, 0.3f, 1.f}, "  Completado.");
            if (!st.result_loaded) {
                ImGui::SameLine();
                if (ImGui::SmallButton("Abrir")) {
                    loadIntoViewer(app, rend, st.result_csv);
                    st.result_loaded = true;
                    app.active_tab = 1;
                }
            }
        } else {
            ImGui::TextColored({1.f, 0.3f, 0.3f, 1.f}, "  Error — ver log.");
        }
    }

    // ── Log ───────────────────────────────────────────────────────────
    ImGui::Separator();
    ImGui::TextDisabled("LOG");
    float avail = ImGui::GetContentRegionAvail().y - 4;
    ImGui::BeginChild("##log", {-1, avail > 40 ? avail : 40}, false,
                      ImGuiWindowFlags_HorizontalScrollbar);
    {
        std::lock_guard<std::mutex> lk(st.log_mutex);
        for (const auto& line : st.log) {
            bool is_err = line.find("[ERROR]") != std::string::npos
                       || line.find("Fatal")   != std::string::npos;
            bool is_ok  = line.find("Done")    != std::string::npos
                       || line.find("\xe2\x86\x92") != std::string::npos; // "→"
            if      (is_err) ImGui::TextColored({1.f,0.4f,0.4f,1.f}, "%s", line.c_str());
            else if (is_ok)  ImGui::TextColored({0.4f,1.f,0.4f,1.f}, "%s", line.c_str());
            else             ImGui::TextUnformatted(line.c_str());
        }
        if (st.scroll_log) { ImGui::SetScrollHereY(1.f); st.scroll_log = false; }
    }
    ImGui::EndChild();
}

// ── Visualization tab ─────────────────────────────────────────────────────────
static void drawVisualizationUI(App& app, Renderer& rend) {
    if (!app.filename.empty()) {
        ImGui::TextDisabled("%s", app.filename.c_str());
        ImGui::Text("%d átomos  /  %d visibles",
                    (int)app.all_atoms.size(), (int)app.visible_indices.size());
    } else {
        ImGui::TextColored({1.f,0.6f,0.3f,1.f}, "Sin datos.");
        ImGui::TextWrapped("Usá la pestaña Analizar o pasá un CSV como argumento.");
        return;
    }
    ImGui::Separator();

    // ── Distortion filter ─────────────────────────────────────────────────
    ImGui::TextDisabled("DISTORTION FILTER");
    bool changed = false;
    changed |= ImGui::SliderFloat("lo", &app.filter_lo, app.dist_min, app.dist_max, "%.4f");
    changed |= ImGui::SliderFloat("hi", &app.filter_hi, app.dist_min, app.dist_max, "%.4f");
    if (app.filter_lo > app.filter_hi) std::swap(app.filter_lo, app.filter_hi);

    // Percentile quick-set buttons
    ImGui::TextDisabled("quick filter:");
    ImGui::SameLine();
    if (ImGui::SmallButton("top 10%%")) {
        float cut = app.dist_min + (app.dist_max - app.dist_min) * 0.9f;
        app.filter_lo = cut; app.filter_hi = app.dist_max; changed = true;
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("top 25%%")) {
        float cut = app.dist_min + (app.dist_max - app.dist_min) * 0.75f;
        app.filter_lo = cut; app.filter_hi = app.dist_max; changed = true;
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("all")) {
        app.filter_lo = app.dist_min; app.filter_hi = app.dist_max; changed = true;
    }
    ImGui::Separator();

    // ── Defect type filter ────────────────────────────────────────────────
    ImGui::TextDisabled("DEFECT TYPE");
    const char* dnames[] = {"Lattice","Interstitial","VacancyAdj","TypeA","Unknown"};
    float cols[5][3] = {
        {0.45f,0.75f,0.45f}, {0.95f,0.30f,0.25f}, {0.25f,0.55f,0.95f},
        {0.95f,0.75f,0.15f}, {0.70f,0.35f,0.90f}
    };
    for (int i = 0; i < 5; ++i) {
        ImGui::PushStyleColor(ImGuiCol_CheckMark, ImVec4(cols[i][0],cols[i][1],cols[i][2],1));
        changed |= ImGui::Checkbox(dnames[i], &app.defect_filter[i]);
        ImGui::PopStyleColor();
    }
    ImGui::Separator();

    // ── Atom type filter ──────────────────────────────────────────────────
    if (!app.atom_types.empty()) {
        ImGui::TextDisabled("ATOM TYPE");
        for (int t : app.atom_types) {
            std::string label = "Type " + std::to_string(t);
            changed |= ImGui::Checkbox(label.c_str(), &app.type_visible[t]);
        }
        ImGui::Separator();
    }

    // ── Color mode ────────────────────────────────────────────────────────
    ImGui::TextDisabled("COLOR BY");
    const char* cmodes[] = {"dist_to_ref","defect_type","atom_type","defect_prob"};
    int cm = (int)app.color_mode;
    if (ImGui::Combo("##color", &cm, cmodes, 4)) {
        app.color_mode = (ColorMode)cm;
        changed = true;
    }
    ImGui::Separator();

    // ── Appearance ────────────────────────────────────────────────────────
    ImGui::TextDisabled("APPEARANCE");
    changed |= ImGui::SliderFloat("radius",  &app.atom_radius, 0.2f, 5.0f, "%.2f Å");
    changed |= ImGui::SliderFloat("opacity", &app.opacity,     0.1f, 1.0f, "%.2f");
    ImGui::Separator();

    // ── Selected atom info ────────────────────────────────────────────────
    if (app.selected_idx >= 0 && app.selected_idx < (int)app.visible_indices.size()) {
        const auto& a = app.all_atoms[app.visible_indices[app.selected_idx]];
        ImGui::TextDisabled("SELECTED ATOM");
        ImGui::Text("id     : %d",    a.id);
        ImGui::Text("type   : %d",    a.type);
        ImGui::Text("pos    : %.3f  %.3f  %.3f", a.x, a.y, a.z);
        ImGui::Text("dist   : %.6f", a.dist_to_ref);
        ImGui::Text("prob   : %.6f", a.defect_prob);
        ImGui::Text("defect : %s",   defectName(a.defect_label));
    }

    if (changed) buildGpuData(app, rend);
}

// ── Main UI window with tab bar ───────────────────────────────────────────────
static void drawUI(App& app, Renderer& rend) {
    ImGui::SetNextWindowPos({0, 0}, ImGuiCond_Always);
    ImGui::SetNextWindowSize({340, (float)app.win_h}, ImGuiCond_Always);
    ImGui::Begin("DistributionTool",
                 nullptr,
                 ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize |
                 ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoTitleBar);

    // Tab bar
    ImGuiTabBarFlags tab_flags = ImGuiTabBarFlags_None;
    if (ImGui::BeginTabBar("##tabs", tab_flags)) {

        // Force tab selection from code (e.g. after auto-load)
        ImGuiTabItemFlags ana_flags  = (app.active_tab == 0) ? ImGuiTabItemFlags_SetSelected : 0;
        ImGuiTabItemFlags vis_flags  = (app.active_tab == 1) ? ImGuiTabItemFlags_SetSelected : 0;
        app.active_tab = -1;  // clear after applying once

        if (ImGui::BeginTabItem("Analizar", nullptr, ana_flags)) {
            ImGui::Spacing();
            drawAnalysisUI(app, rend);
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Visualizar", nullptr, vis_flags)) {
            ImGui::Spacing();
            drawVisualizationUI(app, rend);
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }

    ImGui::End();
}

// ═══════════════════════════════════════════════════════════════════════════════
// main
// ═══════════════════════════════════════════════════════════════════════════════

int main(int argc, char* argv[]) {
    if (!glfwInit()) { std::cerr << "glfwInit failed\n"; return 1; }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    App app;
    Renderer rend;
    g_app  = &app;
    g_rend = &rend;

    GLFWwindow* win = glfwCreateWindow(app.win_w, app.win_h,
                                       "DistributionTool Viewer", nullptr, nullptr);
    if (!win) { std::cerr << "glfwCreateWindow failed\n"; glfwTerminate(); return 1; }
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "gladLoadGL failed\n"; return 1;
    }
    glEnable(GL_MULTISAMPLE);  // activate MSAA requested via GLFW_SAMPLES

    glfwSetScrollCallback(win,          cbScroll);
    glfwSetMouseButtonCallback(win,     cbMouseBtn);
    glfwSetCursorPosCallback(win,       cbCursorPos);
    glfwSetFramebufferSizeCallback(win, cbFramebuffer);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(win, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    // Dark theme tweaks
    auto& style = ImGui::GetStyle();
    style.WindowRounding = 0;
    style.FrameRounding  = 3;
    auto* c = style.Colors;
    c[ImGuiCol_WindowBg]  = {0.10f, 0.10f, 0.12f, 1.f};
    c[ImGuiCol_FrameBg]   = {0.18f, 0.18f, 0.22f, 1.f};
    c[ImGuiCol_SliderGrab]= {0.35f, 0.55f, 0.85f, 1.f};

    rend.init();
    {
        int fbw, fbh;
        glfwGetFramebufferSize(win, &fbw, &fbh);
        rend.initPickFBO(fbw, fbh);
    }

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // ── Load file from argv ───────────────────────────────────────────────────
    if (argc >= 2) {
        try {
            loadIntoViewer(app, rend, argv[1]);
            app.active_tab = 1;  // open straight on Visualizar
            std::cout << "Loaded " << app.all_atoms.size()
                      << " atoms from " << argv[1] << "\n";
        } catch (const std::exception& e) {
            std::cerr << "Load error: " << e.what() << "\n";
        }
    }

    // ── Render loop ───────────────────────────────────────────────────────────
    while (!glfwWindowShouldClose(win)) {
        glfwPollEvents();

        int fbw, fbh;
        glfwGetFramebufferSize(win, &fbw, &fbh);
        app.win_w = fbw; app.win_h = fbh;

        glViewport(0, 0, fbw, fbh);
        glClearColor(0.07f, 0.07f, 0.09f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        float aspect = fbw > 0 ? (float)fbw / (float)fbh : 1.f;
        auto mv   = app.camera.view();
        auto proj = app.camera.proj(aspect);

        rend.draw(mv, proj);

        // ── ImGui ─────────────────────────────────────────────────────────────
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        drawUI(app, rend);

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        // Auto-load result when analysis finishes
        {
            auto& st = app.analysis;
            if (st.done && st.success && st.auto_load && !st.result_loaded) {
                try {
                    loadIntoViewer(app, rend, st.result_csv);
                    st.result_loaded = true;
                    app.active_tab = 1;  // switch to Visualizar
                } catch (const std::exception& e) {
                    std::cerr << "Auto-load failed: " << e.what() << "\n";
                }
            }
        }

        glfwSwapBuffers(win);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(win);
    glfwTerminate();
    return 0;
}
