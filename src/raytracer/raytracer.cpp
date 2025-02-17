#include "raytracer.h"
#include "raytracescene.h"
#include "utils/imagereader.h"
#include <cmath>


RayTracer::RayTracer(Config config) :
    m_config(config)
{}


RGBA toRGBA(const glm::vec4 &illumination) {
    // Task 1
    glm::vec4 clampedIllum = glm::clamp(illumination, 0.0f, 1.0f);

    return RGBA{static_cast<uint8_t>(clampedIllum[0] * 255.0f), static_cast<uint8_t>(clampedIllum[1] * 255.0f), static_cast<uint8_t>(clampedIllum[2] * 255.0f)};
}

void getSphereUV(const glm::vec3 &hitPoint, float &u, float &v) {
    float radius = 0.5f;

    float phi = asin(hitPoint[1] / radius);
    v = (phi) / M_PI + 0.5f;

    float theta = atan2(hitPoint.z,  hitPoint.x);

    if (theta < 0) {
        u = -theta / 2.0f / M_PI;
    }
    else {
        u = 1 - (theta / 2.0f / M_PI);
    }

    if (v == 0 || v == 1) {
        u = 0.5f;
    }

    u = glm::clamp(u, 0.0f, 1.0f);
    v = glm::clamp(v, 0.0f, 1.0f);
    float temp_u = v;
    float temp_v = u;
    u = temp_v;
    v = temp_u;
}

void getCylinderorConeUV(const glm::vec3 &hitPoint, float &u, float &v) {

    v = hitPoint.y + 0.5f;

    float theta = atan2(hitPoint.z,  hitPoint.x);

    if (theta < 0) {
        u = -theta / 2.0f / M_PI;
    }
    else {
        u = 1 - (theta / 2.0f / M_PI);
    }

    u = glm::clamp(u, 0.0f, 1.0f);
    v = glm::clamp(v, 0.0f, 1.0f);

}

void getCubeUV(const glm::vec3 &hitPoint, float &u, float &v) {
    float halfSize = 0.5f;

    if (fabs(hitPoint.x - halfSize) < 1e-3) {
        u = (-hitPoint.z + halfSize) / 1.0f;
        v = (hitPoint.y + halfSize) / 1.0f;
    } else if (fabs(hitPoint.x + halfSize) < 1e-3) {
        u = (hitPoint.z + halfSize) / 1.0f;
        v = (hitPoint.y + halfSize) / 1.0f;
    } else if (fabs(hitPoint.y - halfSize) < 1e-3) {
        u = (hitPoint.x + halfSize) / 1.0f;
        v = (-hitPoint.z + halfSize) / 1.0f;
    } else if (fabs(hitPoint.y + halfSize) < 1e-3) {
        u = (hitPoint.x + halfSize) / 1.0f;
        v = (hitPoint.z + halfSize) / 1.0f;
    } else if (fabs(hitPoint.z - halfSize) < 1e-3) {
        u = (hitPoint.x + halfSize) / 1.0f;
        v = (hitPoint.y + halfSize) / 1.0f;
    } else if (fabs(hitPoint.z + halfSize) < 1e-3) {
        u = (-hitPoint.x + halfSize) / 1.0f;
        v = (hitPoint.y + halfSize) / 1.0f;
    }

    u = glm::clamp(u, 0.0f, 1.0f);
    v = glm::clamp(v, 0.0f, 1.0f);
}

bool intersect (const glm::vec4 ray_origin, const glm::vec4 ray_dir, const RenderShapeData &shape, float &t, glm::vec3 &normal, float &u, float &v) {
    glm::vec4 obj_dir = glm::inverse(shape.ctm) * ray_dir;
    glm::vec4 obj_origin = glm::inverse(shape.ctm) * ray_origin;
    glm::vec3 o = glm::vec3(obj_origin);
    glm::vec3 d = glm::vec3(obj_dir);
    if (shape.primitive.type == PrimitiveType::PRIMITIVE_SPHERE) {
        float radius = 0.5f;
        float a = glm::dot(d, d);
        float b = 2.0f * glm::dot(o, d);
        float c = glm::dot(o, o) - radius * radius;
        float discriminant = b * b - 4 * a * c;
        bool hit;

        if (discriminant < 0) {
            hit = false;
        }
        else {
            float t1 = (-b - sqrt(discriminant)) / (2.0f * a);
            float t2 = (-b + sqrt(discriminant)) / (2.0f * a);

            if (t1 < 0 && t2 < 0) {
                hit = false;
            }
            else if (t1 > 0 && t2 < 0) {
                t = t1;
                hit = true;
            }
            else if (t1 < 0 && t2 > 0) {
                hit = true;
                t = t2;
            }
            else {
                if (t1 <= t2) {
                    hit = true;
                    t = t1;
                }
                else {
                    hit = true;
                    t = t2;
                }
            }
            glm::vec3 hitPoint = o + t * d;
            getSphereUV(hitPoint, u, v);
            normal = glm::normalize(2.0f * hitPoint);
        }
        return hit;
    }
    if (shape.primitive.type == PrimitiveType::PRIMITIVE_CYLINDER) {
        bool top = false;
        float radius = 0.5f;
        bool hit = false;

        float a = d.x * d.x + d.z * d.z;
        float b = 2.0f * (o.x * d.x + o.z * d.z);
        float c = o.x * o.x + o.z * o.z - radius * radius;

        float discriminant = b * b - 4.0f * a * c;

        float tSide = INFINITY;
        float tCaps = INFINITY;
        glm::vec3 normalSide;
        glm::vec3 normalCaps;

        if (discriminant >= 0.0f) {
            float t0 = (-b - sqrt(discriminant)) / (2.0f * a);
            float t1 = (-b + sqrt(discriminant)) / (2.0f * a);

            for (float tCandidate : {t0, t1}) {
                if (tCandidate > 0.0f && tCandidate < t) {
                    float y = o.y + tCandidate * d.y;
                    if (y >= -0.5f && y <= 0.5f) {
                        tSide = tCandidate;
                        glm::vec3 hitPoint = o + tSide * d;
                        normalSide = glm::normalize(glm::vec3(hitPoint.x, 0.0f, hitPoint.z));
                        hit = true;
                        break;
                    }
                }
            }
        }
        if (d.y != 0.0f) {
            float tTop = (0.5f - o.y) / d.y;
            if (tTop > 0.0f && tTop < t) {
                glm::vec3 hitPoint = o + tTop * d;
                float dist2 = hitPoint.x * hitPoint.x + hitPoint.z * hitPoint.z;
                if (dist2 <= radius * radius) {
                    tCaps = tTop;
                    normalCaps = glm::vec3(0.0f, 1.0f, 0.0f);
                    hit = true;
                    top = true;
                }
            }
            float tBot = (-0.5f - o.y) / d.y;
            if (tBot > 0.0f && tBot < t) {
                glm::vec3 hitPoint = o + tBot * d;
                float dist2 = hitPoint.x * hitPoint.x + hitPoint.z * hitPoint.z;
                if (dist2 <= radius * radius && tBot < tCaps) {
                    tCaps = tBot;
                    normalCaps = glm::vec3(0.0f, -1.0f, 0.0f);
                    hit = true;
                    top = false;
                }
            }
        }
        if (tSide < tCaps) {
            if (tSide < t) {
                t = tSide;
                normal = normalSide;
                glm::vec3 hitpoint = o + t * d;
                getCylinderorConeUV(hitpoint, u, v);

                u = glm::clamp(u, 0.0f, 1.0f);
                v = glm::clamp(v, 0.0f, 1.0f);
            } else {
                hit =  false;
            }
        } else {
            if (tCaps < t) {
                t = tCaps;
                normal = normalCaps;
                glm::vec3 hitpoint = o + t * d;
                u = hitpoint[0] + 0.5f;
                if (top == false) {

                    v = hitpoint[2] + 0.5f;
                }
                else {
                    v = -hitpoint[2] + 0.5f;
                }

                u = glm::clamp(u, 0.0f, 1.0f);
                v = glm::clamp(v, 0.0f, 1.0f);
            } else {
                hit = false;
            }
        }

        return hit;

    }
    if (shape.primitive.type == PrimitiveType::PRIMITIVE_CUBE) {
        float tMin = -INFINITY;
        float tMax = INFINITY;
        glm::vec3 tMinNormal;
        glm::vec3 tMaxNormal;

        for (int i = 0; i < 3; ++i) {
            float origin = o[i];
            float direction = d[i];
            float t1, t2;
            glm::vec3 normal1(0.0f);
            glm::vec3 normal2(0.0f);

            if (direction != 0.0f) {
                t1 = (-0.5f - origin) / direction;
                t2 = (0.5f - origin) / direction;

                normal1[i] = -1.0f;
                normal2[i] = 1.0f;
            } else {
                if (origin < -0.5f || origin > 0.5f) {
                    return false;
                }
                t1 = -INFINITY;
                t2 = INFINITY;
            }

            if (t1 > t2) {
                std::swap(t1, t2);
                std::swap(normal1, normal2);
            }
            if (t1 > tMin) {
                tMin = t1;
                tMinNormal = normal1;
            }
            if (t2 < tMax) {
                tMax = t2;
                tMaxNormal = normal2;
            }

            if (tMin > tMax) {
                return false;
            }
            if (tMax < 0.0f) {
                return false;
            }
        }

        if (tMin > 0.0f && tMin < t) {
            t = tMin;
            normal = tMinNormal;
            glm::vec3 hitpoint = o + t * d;
            getCubeUV(hitpoint, u, v);
        } else if (tMax > 0.0f && tMax < t) {
            t = tMax;
            normal = tMaxNormal;
            glm::vec3 hitpoint = o + t * d;
            getCubeUV(hitpoint, u, v);
        } else {
            return false;
        }

        return true;


    }
    if (shape.primitive.type == PrimitiveType::PRIMITIVE_CONE) {
        float height = 1.0f;
        float radius = 0.5f;
        float y_apex = 0.5f;
        float y_base = -0.5f;
        bool hit = false;

        float k = radius / height;
        float k2 = k * k;

        float A = d.x * d.x + d.z * d.z - k2 * d.y * d.y;
        float B = 2.0f * (o.x * d.x + o.z * d.z - k2 * (o.y - y_apex) * d.y);
        float C = o.x * o.x + o.z * o.z - k2 * (o.y - y_apex) * (o.y - y_apex);


        float discriminant = B * B - 4.0f * A * C;

        float tCone = std::numeric_limits<float>::infinity();
        float tBase = std::numeric_limits<float>::infinity();
        glm::vec3 normalCone;
        glm::vec3 normalBase;

        if (discriminant >= 0.0f) {
            float sqrtDiscriminant = sqrt(discriminant);
            float t0 = (-B - sqrtDiscriminant) / (2.0f * A);
            float t1 = (-B + sqrtDiscriminant) / (2.0f * A);

            for (float tCandidate : {t0, t1}) {
                if (tCandidate > 0.0f && tCandidate < t) {
                    float y = o.y + tCandidate * d.y;
                    if (y >= y_base && y <= y_apex) {
                        tCone = tCandidate;
                        glm::vec3 hitPoint = o + tCone * d;
                        float nx = hitPoint.x;
                        float ny = -k2 * (hitPoint.y - y_apex);
                        float nz = hitPoint.z;
                        normalCone = glm::normalize(glm::vec3(nx, ny, nz));
                        hit = true;
                        break;
                    }
                }
            }
        }

        if (fabs(d.y) > 1e-6) {
            float tCap = (y_base - o.y) / d.y;
            if (tCap > 0.0f && tCap < t) {
                glm::vec3 hitPoint = o + tCap * d;
                float dist2 = hitPoint.x * hitPoint.x + hitPoint.z * hitPoint.z;
                float radiusAtBase = radius;
                if (dist2 <= radiusAtBase * radiusAtBase) {
                    tBase = tCap;
                    normalBase = glm::vec3(0.0f, -1.0f, 0.0f);
                    hit = true;
                }
            }
        }

        if (tCone < tBase) {
            if (tCone < t) {
                t = tCone;
                normal = normalCone;
                glm::vec3 hitpoint = o + t * d;
                getCylinderorConeUV(hitpoint, u, v);
            } else {
                return false;
            }
        } else {
            if (tBase < t) {
                t = tBase;
                normal = normalBase;
                glm::vec3 hitpoint = o + t * d;
                u = hitpoint[0] + 0.5f;
                v = hitpoint[2] + 0.5f;

                u = glm::clamp(u, 0.0f, 1.0f);
                v = glm::clamp(v, 0.0f, 1.0f);
            } else {
                return false;
            }
        }

        return hit;

    }
    return false;
}

glm::vec4 raytrace(const glm::vec4 &origin, const glm::vec4 &direction, const RayTraceScene &scene, int depth) {
    if (depth == 0) {
        return glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
    }

    float t = INFINITY;
    RenderShapeData closestShape;
    float U = 0.0f;
    float V = 0.0f;
    glm::vec3 normal;
    glm::vec3 current_normal;
    bool hit = false;

    // Check Intersection
    for (const RenderShapeData &shape : scene.getsceneShapes()) {
        float current_t = INFINITY;
        float current_u = 0.0f;
        float current_v = 0.0f;
        hit = intersect(origin, direction, shape, current_t, current_normal, current_u, current_v);
        if (hit && current_t > 0 && current_t < t) {
            U = current_u;
            V = current_v;
            t = current_t;
            closestShape = shape;
            normal = current_normal;
        }
    }

    if (t < INFINITY) {

        glm::vec4 color = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);

        color += scene.getGlobalData().ka * closestShape.primitive.material.cAmbient;

        RGBA t_color;
        // Texture
        float b = 0;

        glm::vec3 hitPoint_local = glm::inverse(closestShape.ctm) * origin + t * glm::inverse(closestShape.ctm) * direction;

        if (closestShape.primitive.material.textureMap.isUsed) {
            b = closestShape.primitive.material.blend;

            float m = closestShape.primitive.material.textureMap.repeatU;
            float n = closestShape.primitive.material.textureMap.repeatV;

            int w = closestShape.texture->width;
            int h = closestShape.texture->height;

            int c = int(floor(U * int(m) * w)) % w;
            int r = int(floor((1.0f - V) * int(n) * h))% h;

            if (U == 1) c = w - 1;
            if (V == 0) r = h - 1;

            t_color = closestShape.texture->data[r * w + c];
        }
        else {
            t_color = RGBA{0, 0, 0, 255};
        }
        glm::vec4 texture_color = glm::vec4(t_color.r / 255.0f, t_color.g / 255.0f, t_color.b / 255.0f, 1.0f);

        //

        glm::vec3 hitPoint = glm::vec3(closestShape.ctm * glm::vec4(hitPoint_local, 1.0f));

        normal = glm::normalize(glm::inverse(glm::transpose(glm::mat3(closestShape.ctm))) * normal);


        // Shadow
        for (const SceneLightData &light : scene.getsceneLights()) {
            float shadow = 1.0f;
            float U_notuse;
            float V_notuse;
            float shadow_t = INFINITY;
            RenderShapeData shadowShape;
            glm::vec3 normal_No_use;
            glm::vec3 current_normal_No_use;
            bool shadow_hit = false;
            glm::vec4 shadow_ray_dir;

            if (light.type == LightType::LIGHT_DIRECTIONAL) {
                shadow_ray_dir =  glm::normalize(-light.dir);
            }
            else {
                shadow_ray_dir = glm::normalize(glm::vec4(glm::vec3(light.pos) - hitPoint, 0.0f));
            }
            glm::vec4 shadow_ray_origin = glm::vec4(hitPoint, 1.0f) + 0.01f * shadow_ray_dir;


            for (const RenderShapeData &shape : scene.getsceneShapes()) {
                float current_shadow_t = INFINITY;
                shadow_hit = intersect(shadow_ray_origin, shadow_ray_dir, shape, current_shadow_t, current_normal_No_use, U_notuse, V_notuse);
                if (shadow_hit && current_shadow_t > 0 && current_shadow_t < t) {
                    shadow_t = current_shadow_t;
                    shadowShape = shape;
                }
            }

            if (shadow_t < INFINITY) {
                if (light.type == LightType::LIGHT_DIRECTIONAL) {
                    shadow = 0.0f;
                }
                else {
                    float light_distance = glm::length(glm::vec3(light.pos) - hitPoint);

                    glm::vec4 shadow_hitPoint_local = glm::inverse(shadowShape.ctm) * shadow_ray_origin +
                                                          shadow_t * glm::inverse(shadowShape.ctm) * shadow_ray_dir;
                    glm::vec3 shadow_hitPoint = glm::vec3(shadowShape.ctm * shadow_hitPoint_local);

                    float shadow_distance = glm::length(shadow_hitPoint - hitPoint);

                    if (shadow_distance < light_distance) {
                        shadow = 0.0f;
                    }
                }
            }
            //

            // Light
            glm::vec3 light_dir = glm::normalize(hitPoint - glm::vec3(light.pos));
            float distance = glm::length(glm::vec3(light.pos) - hitPoint);
            float attenuation = glm::min(1.0f, 1.0f / (light.function.x + light.function.y * distance + light.function.z * distance * distance));

            float diffuseFactor = glm::max(glm::dot(normal, -light_dir), 0.0f);

            glm::vec3 reflectDir = glm::reflect(light_dir, normal);
            glm::vec3 dirTocam = glm::normalize(glm::vec3(scene.getCamera().getcameraData().pos) - hitPoint);
            float specularFactor = glm::pow(glm::max(glm::dot(dirTocam, reflectDir), 0.0f), closestShape.primitive.material.shininess);

            if (light.type == LightType::LIGHT_POINT) {


                color[0] += shadow * ((1 - b) * scene.getGlobalData().kd * closestShape.primitive.material.cDiffuse[0] + b * texture_color[0]) * diffuseFactor * attenuation * light.color.r;
                color[1] += shadow * ((1 - b) * scene.getGlobalData().kd * closestShape.primitive.material.cDiffuse[1] + b * texture_color[1]) * diffuseFactor * attenuation * light.color.g;
                color[2] += shadow * ((1 - b) * scene.getGlobalData().kd * closestShape.primitive.material.cDiffuse[2] + b * texture_color[2]) * diffuseFactor * attenuation * light.color.b;

                color[0] += shadow * scene.getGlobalData().ks * closestShape.primitive.material.cSpecular[0] * attenuation * light.color.r * specularFactor;
                color[1] += shadow * scene.getGlobalData().ks * closestShape.primitive.material.cSpecular[1] * attenuation * light.color.g * specularFactor;
                color[2] += shadow * scene.getGlobalData().ks * closestShape.primitive.material.cSpecular[2] * attenuation * light.color.b * specularFactor;
            }

            else if (light.type == LightType::LIGHT_SPOT) {

                float x = acos(glm::dot(light_dir, glm::vec3(glm::normalize(light.dir))));

                float inner = light.angle - light.penumbra;

                float base = (x - inner)/(light.angle - inner);

                float falloff = -2.0f * pow(base, 3.0) + 3.0f * pow(base, 2.0);

                if (x <= inner) {
                    color[0] += shadow * ((1 - b) * scene.getGlobalData().kd * closestShape.primitive.material.cDiffuse[0] + b * texture_color[0]) * diffuseFactor * attenuation * light.color.r;
                    color[1] += shadow * ((1 - b) * scene.getGlobalData().kd * closestShape.primitive.material.cDiffuse[1] + b * texture_color[1]) * diffuseFactor * attenuation * light.color.g;
                    color[2] += shadow * ((1 - b) * scene.getGlobalData().kd * closestShape.primitive.material.cDiffuse[2] + b * texture_color[2]) * diffuseFactor * attenuation * light.color.b;

                    color[0] += shadow * scene.getGlobalData().ks * closestShape.primitive.material.cSpecular[0] * attenuation * light.color.r * specularFactor;
                    color[1] += shadow * scene.getGlobalData().ks * closestShape.primitive.material.cSpecular[1] * attenuation * light.color.g * specularFactor;
                    color[2] += shadow * scene.getGlobalData().ks * closestShape.primitive.material.cSpecular[2] * attenuation * light.color.b * specularFactor;
                }
                else if (x > inner && x <= light.angle){
                    color[0] += shadow * ((1 - b) * scene.getGlobalData().kd * closestShape.primitive.material.cDiffuse[0] + b * texture_color[0]) * diffuseFactor * attenuation * (light.color.r * (1 - falloff));
                    color[1] += shadow * ((1 - b) * scene.getGlobalData().kd * closestShape.primitive.material.cDiffuse[1] + b * texture_color[1]) * diffuseFactor * attenuation * (light.color.g * (1 - falloff));
                    color[2] += shadow * ((1 - b) * scene.getGlobalData().kd * closestShape.primitive.material.cDiffuse[2] + b * texture_color[2]) * diffuseFactor * attenuation * (light.color.b * (1 - falloff));

                    color[0] += shadow * scene.getGlobalData().ks * closestShape.primitive.material.cSpecular[0] * attenuation * (light.color.r * (1 - falloff)) * specularFactor;
                    color[1] += shadow * scene.getGlobalData().ks * closestShape.primitive.material.cSpecular[1] * attenuation * (light.color.g * (1 - falloff)) * specularFactor;
                    color[2] += shadow * scene.getGlobalData().ks * closestShape.primitive.material.cSpecular[2] * attenuation * (light.color.b * (1 - falloff)) * specularFactor;
                }
            }
            else if (light.type == LightType::LIGHT_DIRECTIONAL) {
                float diffuseFactor = glm::max(glm::dot(normal, glm::vec3(-light.dir)), 0.0f);

                color[0] += shadow * ((1 - b) * scene.getGlobalData().kd * closestShape.primitive.material.cDiffuse[0] + b * texture_color[0]) * diffuseFactor * light.color.r;
                color[1] += shadow * ((1 - b) * scene.getGlobalData().kd * closestShape.primitive.material.cDiffuse[1] + b * texture_color[1]) * diffuseFactor * light.color.g;
                color[2] += shadow * ((1 - b) * scene.getGlobalData().kd * closestShape.primitive.material.cDiffuse[2] + b * texture_color[2]) * diffuseFactor * light.color.b;

                glm::vec3 reflectDir = glm::reflect(glm::vec3(light.dir), normal);
                glm::vec3 dirTocam = glm::normalize(glm::vec3(scene.getCamera().getcameraData().pos) - hitPoint);
                float specularFactor = glm::pow(glm::max(glm::dot(dirTocam, reflectDir), 0.0f), closestShape.primitive.material.shininess);
                color[0] += shadow * scene.getGlobalData().ks * closestShape.primitive.material.cSpecular[0] * light.color.r * specularFactor;
                color[1] += shadow * scene.getGlobalData().ks * closestShape.primitive.material.cSpecular[1] * light.color.g * specularFactor;
                color[2] += shadow * scene.getGlobalData().ks * closestShape.primitive.material.cSpecular[2] * light.color.b * specularFactor;
            }
        }
        //

        // Reflection
        if (closestShape.primitive.material.cReflective != glm::vec4(0)) {
            glm::vec3 directionToOrigin = glm::normalize(glm::vec3(origin) - hitPoint);
            glm::vec4 reflected_dir = glm::vec4((2.0f * glm::dot(normal, directionToOrigin) * normal - directionToOrigin), 0.0f);
            glm::vec4 reflected_origin = glm::vec4(hitPoint, 1.0f) + 0.01f * reflected_dir;

            glm::vec4 reflected_color = raytrace(reflected_origin, reflected_dir, scene, depth - 1);
            color[0] += scene.getGlobalData().ks * closestShape.primitive.material.cReflective[0] * reflected_color[0];
            color[1] += scene.getGlobalData().ks * closestShape.primitive.material.cReflective[1] * reflected_color[1];
            color[2] += scene.getGlobalData().ks * closestShape.primitive.material.cReflective[2] * reflected_color[2];

        }
        //
        return color;
    }

    return glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
}

void RayTracer::render(RGBA *imageData, const RayTraceScene &scene) {
    int height = scene.height();
    int width = scene.width();

    float aspectRatio = static_cast<float>(width) / static_cast<float>(height);
    float vh = 2 * tan(scene.getCamera().getHeightAngle() / 2.0);
    float vw = vh * aspectRatio;

    for (int i = 0; i < height; i++){
        for (int j = 0; j < width; j++) {
            float y = vh * (((height - 1 - i + 0.5) / height) - 0.5);
            float x = vw * (((j + 0.5) / width) - 0.5);
            float z = -1.0;

            glm::vec4 pos = glm::vec4(x, y ,z, 1);
            glm::vec4 eye = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
            glm::vec4 d = glm::normalize(scene.getCamera().getInverseViewMatrix() * glm::normalize(pos - eye));
            glm::vec4 origin = scene.getCamera().getInverseViewMatrix() * eye;

            glm::vec4 color = raytrace(origin, d, scene, 4);
            imageData[i * width + j] = toRGBA(color);
        }
    }
}
