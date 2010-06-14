sampler2D myTexture;
varying vec2 vTexCoord;

void main (void)  
{  vec4 color = texture2D(myTexture, vTexCoord);       
   
   if (color.rgb == vec3(1.0,0.0,0.0))
      discard; 
   
   gl_FragColor = color;
}    
