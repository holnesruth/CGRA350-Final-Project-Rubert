#version 330 core

#define M_PI 3.1415926535897932384626433832795
#define ETA 1.4

#define RED 700
#define GREEN 560
#define BLUE 460

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform mat4 uModelMatrix;

uniform sampler2D uTexture;
uniform sampler2D uNoise;
uniform samplerCube uCubeMap;

uniform vec2 uThickness;
uniform vec2 uLightEffects;

uniform float uTime;
uniform vec2 uResolution;

uniform bool uFlow;

// viewspace data (this must match the output of the vertex shader)
in VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
} f_in;

// framebuffer output
out vec4 fb_color;

const vec3 DIRECTION = vec3(1, 1, 4);

/*============================ FLOW NOISE IMPLEMENTATION ============================*/

float time = uTime*0.1;

float noise(vec2 x) {
	return texture(uNoise, x*0.1).x; 
}

float flow(vec2 v) {
	float divisor = 2.0;
	float result = 0.0;
	vec2 base = v;
	float octaves = 7.0;

	for (float i = 1.0; i < octaves; i++) {
		//primary flow speed
		v += time*0.6;
		
		//add noise octave
		result += (sin(noise(v)*octaves)*0.5 + 0.5)/divisor;
		
		//blend factor - 0.5 being low, 0.95 being high
		v = mix(base, v, 0.77);
		
		//intensity scaling
		divisor *= 1.4;

		//octave scaling
		v *= 2.0;
		base *= 1.9;
	}

	return result;	
}

/*============================ REFRACTION IMPLEMENTATION ============================*/

float snellsLaw(float theta_i, float n1, float n2) {
	return asin((n1/n2)*sin(theta_i));
}

float fresnel(float theta_i, float theta_t, float n1, float n2) {
	float num = n1*cos(theta_i) - n2*cos(theta_t);
	float den = n1*cos(theta_i) + n2*cos(theta_t);
	float Rs = pow(abs(num/den), 2);

	num = n1*cos(theta_t) - n2*cos(theta_i);
	den = n1*cos(theta_t) + n2*cos(theta_i);
	float Rp = pow(abs(num/den), 2);

	return 2.0*(Rs + Rp); //Scalar should be 0.5 but the colours aren't very vibrant
}

float calculateLightColor(float lambda, float thickness, float nAir, float nFilm, float theta_i, float intensity) {
	float theta_t = max(snellsLaw(theta_i, nAir, nFilm), 0.0);

	float d = 2 * M_PI * thickness * nFilm * cos(theta_t);
	float sind = sin(d/lambda);

	float fresnel = fresnel(theta_i, theta_t, nAir, nFilm);

	return 4.0 * intensity * fresnel * sind * sind;
}

void main() {
	float minThickness = uThickness.x;
	float maxThickness = uThickness.y;
	float intensity = uLightEffects.x;
	float transparency = uLightEffects.y;

	vec3 norm = f_in.normal;

	vec4 light = uProjectionMatrix * vec4(DIRECTION, 1.0);
	vec3 lightDir = normalize(vec3(light));

	vec3 viewDir = normalize(-f_in.position);
	vec3 reflectDir = reflect(-lightDir, norm);

	vec3 thickness = texture(uTexture, f_in.textureCoord).rgb;
	float t = (thickness.x + thickness.y + thickness.z)/3.0;
	float w = ((minThickness*(1.0 - t)) + (maxThickness*t));
	
//	vec2 v = (f_in.position.xy / uResolution.xy) - 0.5;
//	v.x *= uResolution.x/uResolution.y;
//	v *= 3.0;

	//introduce noise
	if (uFlow) {
		w /= flow(f_in.normal.xy); 
	} else {
		w /= noise(f_in.textureCoord.xy);
	}

	float theta_i = max(dot(norm, lightDir), 0.0);

	float red = calculateLightColor(RED, w, 1, ETA, theta_i, intensity);
	float green = calculateLightColor(GREEN, w, 1, ETA, theta_i, intensity);
	float blue = calculateLightColor(BLUE, w, 1, ETA, theta_i, intensity);

	vec3 filmColour = vec3(red, green, blue);
	vec3 reflection = texture(uCubeMap, reflectDir).rgb;

	vec3 finalShading = ((1 - transparency) * filmColour) + (intensity * reflection);
	
	// output to the framebuffer
	fb_color = vec4(finalShading, transparency);
}