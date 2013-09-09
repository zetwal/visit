const int BLEND = 1;                      
 
uniform sampler2D tex7;  // front-to-back  
uniform sampler2D tex8;  // back-to-front  
uniform int option;                        
 
varying vec4 clip_Vertex;                  
 
void main(){ 
    vec4 color, color1, color2;    
    color1 = texture2D(tex7, gl_TexCoord[0].xy);       
 
    if (option == BLEND){          
        color2 = texture2D(tex8, gl_TexCoord[0].xy);   
        color = color1 + color2*(1.0-color1.a);        
    } else{ 
		color = color1;             
    } 
 
	gl_FragColor = color;          
}