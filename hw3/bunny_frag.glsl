#version 330 core

vec3 lightPos = vec3(-5, 5, -5);
vec3 eyePos = vec3(-10, 0, 0);

vec3 I = vec3(0.8, 0.8, 0.1);
vec3 Iamb = vec3(0.8, 0.8, 0.8);

vec3 kd = vec3(0.8, 0.6, 0.3);
// vec3 ka = vec3(0.1, 0.1, 0.1);
vec3 ks = vec3(0.8, 0.6, 0.3);

in vec4 fragPos;
in vec3 N;

out vec4 outColor;
void main(void)
{
	vec3 L = normalize(lightPos - vec3(fragPos));
	vec3 V = normalize(eyePos - vec3(fragPos));
	vec3 H = normalize(L + V);
	float NdotL = dot(N, L);
	float NdotH = dot(N, H);

	vec3 diffuseColor = I * kd * max(0, NdotL);
	// vec3 ambientColor = Iamb * ka;
	vec3 specularColor = I * ks * pow(max(0, NdotH), 20);

    outColor = vec4(diffuseColor + specularColor, 1);
}
