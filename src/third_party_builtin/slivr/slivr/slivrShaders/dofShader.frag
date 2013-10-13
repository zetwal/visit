const int BACK_TO_FRONT = 0;
const int FRONT_TO_BACK = 1;

uniform vec4 loc0, loc1, loc2, loc3, loc4, loc5;
uniform sampler3D tex0; // scalar values - coarse
uniform sampler3D tex1; // gradient - unused

uniform sampler1D tex2;	// 1D transfer function
uniform sampler2D tex3; // 2D transfer function

uniform sampler2D tex6;   //Buffer of last slice   
 
uniform vec2 windowSize;           
uniform int numSamples;            
uniform float angleFactor;                
 
uniform int traversalDirection;    // required for DoF
 
varying vec4 clip_Vertex;          
varying vec2 transformed_c;        

vec4 blur(vec4 source){
    int num_of_samples;

    vec3 samplePosition[9];             
	samplePosition[0] = vec3( 0, 0, 1); 

	samplePosition[1] = vec3( 0, 1, 1);    
   	samplePosition[2] = vec3( 0,-1, 1);    
   	samplePosition[3] = vec3( 1, 0, 1);    
   	samplePosition[4] = vec3(-1, 0, 1);

   	samplePosition[5] = vec3(-1,-1, 1);    
   	samplePosition[6] = vec3(-1, 1, 1);    
   	samplePosition[7] = vec3( 1,-1, 1);    
   	samplePosition[8] = vec3( 1, 1, 1);

    num_of_samples = numSamples;

    // compute blurred part
    vec4 blurred = vec4(0.0);
    for (int i=0; i<num_of_samples; i++){
        vec2 xy_d = (clip_Vertex.xy/clip_Vertex.w);
        vec2 p_d = xy_d + transformed_c * ((angleFactor*samplePosition[i].xy)/(windowSize));
        vec2 p_t = 0.5 * p_d + vec2(0.5);
        
        blurred += texture2D(tex6, p_t);
    }
    
    vec4 destination = blurred/float(num_of_samples);
    return destination;
}



void main()
{
	vec4 l = loc0; 		// {lx, ly, lz, alpha} VOL_LIT_HEAD
	vec4 k = loc1; 		// {ka, kd, ks, ns}
	vec4 g = loc2; 		// {1/gradrange, -gradmin/gradrange, 0, 0}
	vec4 dir = loc4; 	// VOL_GRAD_COMPUTE_2_1

    vec4 t = gl_TexCoord[0];
   	vec4 v = texture3D(tex0, t.stp); 	// VOL_VLUP_1_1
   	vec4 c, n;
   	bool transferFn1D, lightingOn;
	
   	if (loc5.x < 0.5)
   		transferFn1D = true;
   	else
   		transferFn1D = false;
    
   	if (loc5.y < 0.5)
   		lightingOn = false;
   	else
   		lightingOn = true;
    
	// Lighting calculation
	// Compute normal
	if ((lightingOn == true) || (transferFn1D == false)){
		vec4 w, p, r;
		mat4 tmat = gl_TextureMatrixInverseTranspose[0];
        
        v = vec4(v.x); 
		n = vec4(0.0);
		
		w = vec4(0.0);
		w.x = dir.x;
	   	p = clamp(gl_TexCoord[0] + w, 0.0, 1.0);
	   	r = texture3D(tex0, p.stp);
	   	n.x = r.x + n.x;
	   	p = clamp(gl_TexCoord[0] - w, 0.0, 1.0);
	   	r = texture3D(tex0, p.stp);
	   	n.x = r.x - n.x;
        
	   	w = vec4(0.0);
	   	w.y = dir.y;
	   	p = clamp(gl_TexCoord[0] + w, 0.0, 1.0);
	   	r = texture3D(tex0, p.stp);
	   	n.y = r.x + n.y;
	   	p = clamp(gl_TexCoord[0] - w, 0.0, 1.0);
	   	r = texture3D(tex0, p.stp);
	   	n.y = r.x - n.y;
        
		w = vec4(0.0);
	   	w.z = dir.z;
	   	p = clamp(gl_TexCoord[0] + w, 0.0, 1.0);
	   	r = texture3D(tex0, p.stp);
	   	n.z = r.x + n.z;
	   	p = clamp(gl_TexCoord[0] - w, 0.0, 1.0);
	   	r = texture3D(tex0, p.stp);
	   	n.z = r.x - n.z;
        
	   	//      b.recompute the normal
	   	w.x = dot(n.xxx, vec3(tmat[0].x, tmat[1].x, tmat[2].x));
	   	w.y = dot(n.yyy, vec3(tmat[0].y, tmat[1].y, tmat[2].y));
	   	w.z = dot(n.zzz, vec3(tmat[0].z, tmat[1].z, tmat[2].z));
	   	n.xyz = w.xyz;
	   	n = normalize(n);
   	}

	// Query the transfer function
	if (transferFn1D == true)
		c = texture1D(tex2, v.x); 	// VOL_TFLUP_1_1
	else{
		v.y = n.y * 1.75; // 1.75 is unexplained from the original slivr
		c = texture2D(tex3, v.xy); 	// VOL_TFLUP_2_1
	}
    
	// Compute lighting based on normal
	if (lightingOn == true){
		// 		a. angle between normal and light
		n.w = clamp(abs(dot(l.xyz,n.xyz)),0.0,1.0);	// angle between light and normal; clamping & abs: two-sided lighting
        
		// Calculate color using phong shading
		// 			I  * ka	   + I  * kd*abs(cos(angle)) + ks*abs(cos(angle))^ns
		c.xyz  = (c.xyz * k.x) +  (c.xyz * (n.w * k.y))  + k.z * pow(n.w, k.w) * c.w;	// Ia * ks*abs(cos(angle))^ns 
	}

	// get previous value
  	vec4 source, destination, composited;
   	source = c;
    vec2 fclip = (0.5 * clip_Vertex.xy/clip_Vertex.w) + vec2(0.5);

   	//destination = texture2D(tex6, fclip);		// normal
   	destination = blur(texture2D(tex6, fclip));	// dof

   	if (traversalDirection == FRONT_TO_BACK) 
       composited = source * (1.0 - destination.a) + destination;  
   	else    
       composited = destination * (1.0 - source.a) + source;       

   	gl_FragData[0] = composited;   
}
