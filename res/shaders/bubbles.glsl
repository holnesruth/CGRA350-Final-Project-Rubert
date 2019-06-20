#version 330 core

// useful values
#define M_PI 3.1415926535897932384626433832795
#define ETA 1.4

// wavelengths
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

uniform bool uShowNoise;
uniform float uTime;
uniform int uOctaves;
uniform vec2 uFlowSpeeds;

// viewspace data (this must match the output of the vertex shader)
in VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
} f_in;

// framebuffer output
out vec4 fb_color;

// light direction
const vec3 DIRECTION = vec3(1, 1, 4);

/*============================ FLOW NOISE ============================*/

float time = uTime*0.1;


/**
 * Sample a noise texture at a specific uv coordinate
 */
float noise(vec2 uv) {
	return texture(uNoise, 0.1*uv).x; 
}


/**
 * Create a 2x2 rotation matrix from an angle, theta
 */
mat2 rotation_matrix(float theta) {
	float c = cos(theta);
	float s = sin(theta);
	return mat2(c, -s, 
				s, c);
}


/**
 * Find the gradient of a vector using the definition of a derivative
 */
vec2 gradient_vec(vec2 v) {
	float epsilon = 0.09;

	float dvdx = noise(vec2(v.x + epsilon, v.y)) - noise(vec2(v.x - epsilon, v.y));
	float dvdy = noise(vec2(v.x, v.y + epsilon)) - noise(vec2(v.x, v.y - epsilon));

	return vec2(dvdx, dvdy);
}


/**
 * Calculate the flow of the film using a specified number of noise octaves
 */
float flow(vec2 pos, float octaves) {
	float divisor = 2.0;
	float result = 0.0;
	vec2 basePos = pos;

	for (float octave = 1.0; octave < octaves; octave++) {
		
		/* === Calculate flow === */

		pos += time * uFlowSpeeds.x;		//dynamic flow speed 
		basePos += time * uFlowSpeeds.y;	//static flow speed
		
		vec2 gradient = gradient_vec(octave*pos*0.34 + time*0.1);		//find the gradient of the point
		gradient *= rotation_matrix(6.0*time - 2.0*pos.x - 1.2*pos.y);	//rotate the gradient

		pos += gradient*0.5;	//add the rotated gradient to the vector
		
		//add noise for this octave to the result
		result += (sin(noise(pos)*octaves) + 1.0)/(2.0 * divisor);
		
		/* === Next octave === */

		//blend the static and dynamic flows together
		//a blend factor of 0.5 is low, 0.95 is high
		pos = mix(basePos, pos, 0.77);
		
		//make noise less intense for next octave
		divisor *= 1.4;

		//add scaling to the position and base between octaves
		pos *= 2.0;
		basePos *= 1.9;
	}

	return result;
}


/*============================ REFLECTANCE ============================*/


/** 
 * Find the refraction angle of light when traveling between two mediums
 *
 * theta_i: The angle of incidence
 * n1: The index of refraction of the material the light is leaving
 * n2: The index of refraction of the material the light is entering
 */
float snellsLaw(float theta_i, float n1, float n2) {
	return asin((n1/n2)*sin(theta_i));
}


/** 
 * Calculate the reflectance of a specific type of polarised light 
 * when travelling between two mediums.
 * 
 * For s-polarised: theta_1 = theta_i and theta_2 = theta_t
 * For p-polarised: theta_1 = theta_t and theta_2 = theta_i
 *
 * theta_1: The first angle in the ratio
 * theta_2: The second angle in the ratio
 * n1: The index of refraction of the material the light is leaving
 * n2: The index of refraction of the material the light is entering
 */
float polarReflectance(float n1, float theta_1, float n2, float theta_2) {
	float minus = n1*cos(theta_1) - n2*cos(theta_2);
	float add = n1*cos(theta_1) + n2*cos(theta_2);

	return pow(abs(minus/add), 2);
}


/** 
 * Calculate the effective reflectance of light when travelling between two mediums
 *
 * theta_i: The angle of incidence
 * theta_t: The angle of refraction
 * n1: The index of refraction of the material the light is leaving
 * n2: The index of refraction of the material the light is entering
 */
float fresnel(float theta_i, float theta_t, float n1, float n2) {
	// s-polarised light
	float Rs = polarReflectance(n1, theta_i, n2, theta_t);

	// p-polarised light
	float Rp = polarReflectance(n1, theta_t, n2, theta_i);

	// effective reflectance
	return 0.5*(Rs + Rp);
}


/** 
 * Calculate the colour of light that reflects from the surface
 *
 * lambda: The wavelength of the light
 * thickness: The thickness of the film
 * n1: The index of refraction of the material the light is leaving
 * n2: The index of refraction of the material the light is entering
 * theta_i: The angle of incidence
 * intensity: The intensity of the light ray
 */
float calculateLightColor(float lambda, float thickness, float n1, float n2, float theta_i, float intensity) {
	float theta_t = max(snellsLaw(theta_i, n1, n2), 0.0);
	float reflectance = 4.0 * fresnel(theta_i, theta_t, n1, n2);

	// Extra distance covered by rays that travel through the film before reflecting
	float d = 2.0 * M_PI * thickness * n2 * cos(theta_t);
	float sin_d = sin(d/lambda);

	return 4.0 * intensity * reflectance * sin_d * sin_d;
}


/* ========================= MAIN ======================== */


void main() {
	//get the parameters out of the uniforms
	float minThickness = uThickness.x;
	float maxThickness = uThickness.y;
	float intensity = uLightEffects.x;
	float transparency = uLightEffects.y;

	//save the normal vector for later
	vec3 norm = f_in.normal;

	// find the thickness of the film for this fragment
	vec3 thickness = texture(uTexture, f_in.textureCoord).rgb;
	float t = (thickness.x + thickness.y + thickness.z)/3.0;
	float w = (1.0 - t)*minThickness + t*maxThickness;
	if (uShowNoise) w /= flow(norm.xy, uOctaves);	//introduce noise

	//find the angle of incidence of the light
	vec4 light = uProjectionMatrix * vec4(DIRECTION, 1.0);
	vec3 lightDir = normalize(vec3(light));
	float theta_i = max(dot(norm, lightDir), 0.0);

	//calculate the colour of reflected light for rgb wavelengths
	float red = calculateLightColor(RED, w, 1.0, ETA, theta_i, intensity);
	float green = calculateLightColor(GREEN, w, 1.0, ETA, theta_i, intensity);
	float blue = calculateLightColor(BLUE, w, 1.0, ETA, theta_i, intensity);

	vec3 filmColour = vec3(red, green, blue);

	//calculate the reflection for this fragment
	vec3 reflectDir = reflect(-lightDir, norm);
	vec3 reflection = texture(uCubeMap, reflectDir).rgb;

	//combine the film colour and reflections for the result
	vec3 finalShading = (1 - transparency)*filmColour + intensity*reflection;
	
	// output to the framebuffer
	fb_color = vec4(finalShading, transparency);
}