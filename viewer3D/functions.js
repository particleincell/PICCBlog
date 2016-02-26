var results_3d;

function init()
{                 
	log_ele = document.getElementById("log");
	
	results_3d = new Results3D();
	
	/*set up drag and drop*/
	setDropListener(document.getElementById("drop"));  
	results_3d.tick();
	
}
  
/*adds listener for drop events*/
function setDropListener(dropZone)
{
	dropZone.addEventListener('dragover', handleDragOver, false);
	dropZone.addEventListener('drop', handleFileSelect, false);
	dropZone.addEventListener('dragleave',handleDragLeave, true);
}

/*****************************************************/

var log_ele;
var interval_id;

function log(text) {log_ele.innerHTML+=text+"<br>";}
function log_error(text) {log_ele.innerHTML+="<span class='red'>"+text+"</span><br>";}
function log_clear(){log_ele.innerHTML="";}


function isNumber(val) {return isFinite(parseFloat(val));}

/********** 3D ******************************************/

 var CLICK = {
	  LEFT : {value: 0, name: "Left"}, 
	  RIGHT: {value: 1, name: "Right"}, 
	  MIDDLE : {value: 2, name: "Middle"}
	};

function degToRad(d) {return d*Math.PI/180.0;}
	
var Results3D = function()
{
	var p = Results3D.prototype;
	
	p.init = function()
	{
		/*get file name from storage*/
		var file_name = "surface.vtu";
			
		document.getElementById("c_3d").addEventListener("contextmenu", function(evt) {evt.preventDefault();});
		this.canvas = new Canvas3D("c_3d");
		
		/*hooks for variable data*/
		document.getElementById("vars_3d").addEventListener("change", function(t){return function(evt){var index=evt.target.selectedIndex;t.changeVar(index);}}(this));
		var f_range = function(t){return function(e){t.changeRange(e);}}(this);
		document.getElementById("range_min_3d").addEventListener("input", f_range);
		document.getElementById("range_max_3d").addEventListener("input", f_range);
		document.getElementById("range_levels_3d").addEventListener("input", f_range);
		document.getElementById("range_log_3d").addEventListener("change", f_range);
		
		this.components = [{}];	/*first component is for geometry*/
		this.var_index = 0;	/*current variable*/
	}
	
	/*sets data to stuff from file*/
	p.setData = function(data)
	{
		this.canvas.setComponent(data,this.components[0]);
		
		this.initVars();
		this.canvas.convertScalars(this.var_index);
		this.canvas.resetView(true);
		
	}
	
	/*changes current var*/
	p.changeVar = function(index)
	{
		this.var_index=index;

		var comp = this.components[0];
		document.getElementById("range_min_3d").value = comp.scalars[index].min;
		document.getElementById("range_max_3d").value = comp.scalars[index].max;
		document.getElementById("range_levels_3d").value = comp.scalars[index].levels;

		this.canvas.convertScalars(index);
	}

	/*changes limits*/
	p.changeRange = function(evt)
	{
		var comp = this.components[0];
		
		comp.scalars[this.var_index].min = document.getElementById("range_min_3d").value;
		comp.scalars[this.var_index].max = document.getElementById("range_max_3d").value;
		comp.scalars[this.var_index].levels = document.getElementById("range_levels_3d").value;
		comp.scalars[this.var_index].log = document.getElementById("range_log_3d").checked;
		this.canvas.convertScalars(this.var_index);
	}
  
   /*resizes canvas div*/
	p.checkSize = function()
	{
		var DELTA=15;	/*need to add little buffer so the resize handle shows up*/
		var cns=[this.canvas];
		for (var i=0;i<cns.length;i++)
		{
			var canvas = cns[i].canvas;
			var ele = canvas.parentElement;
			var bw = ele.offsetWidth;
			var bh = ele.offsetHeight;
			if ((bw>0 && bh>0) && (bw!=(cns[i].w+DELTA) || bh!=(cns[i].h+DELTA)))
			{
				cns[i].resize(bw-DELTA,bh-DELTA);
			}
		}	
	}

	/*fills select box with loaded field data*/
	p.initVars = function()
	{
		var list = document.getElementById("vars_3d");
		var v=0;
		/*grab currently selected variable*/
		var length=list.length;	/*save length since removing will decrease this*/
		for (var i=0;i<list.length;i++) {if (list[i].selected) v=i;}
		for (var i=0;i<length;i++) {list.remove(0);}	/*remove first*/
		
		var comp = this.components[0];
		for (var i=0;i<comp.scalars.length;i++)
		{
			var name = comp.scalars[i].name;
			var option = document.createElement("option");
			option.text = name;
			option.value = i;
			if (i==v) option.selected=true;
			list.add(option);
		}
		
		this.changeVar(v);
	}
	
	p.tick = function() 
	{
		this.checkSize();
		this.canvas.drawScene();	
		requestAnimationFrame(function(t){return function(){t.tick();}}(this));		
	}

	this.init();	
}

/********************************************************************************************/
/************** CANVAS 3D *******************************************************************/
/*plotting in an associated canvas*/
var Canvas3D = function(id)
{
	var p=Canvas3D.prototype;
	p.init = function(id) 
	{
		this.canvas=document.getElementById(id);
		this.ctx = this.canvas.getContext("webgl");
		this.w = this.canvas.width;
		this.h = this.canvas.height;
		
		/*init webgl*/
		this.initGL();
		this.initShaders()
        
		var f_up=function(t){return function(e){t.handleMouseUp(e);}}(this);
		var f_move=function(t){return function(e){t.handleMouseMove(e);}}(this);
		
		this.canvas.addEventListener("mousedown",function(t){return function(e){t.handleMouseDown(e);}}(this));
		this.canvas.addEventListener("mouseup",f_up);
		this.canvas.addEventListener("mousemove",f_move);
		
		window.addEventListener("mouseup",f_up);
        window.addEventListener("mousemove", f_move);
		
		var gl = this.gl;
		gl.clearColor(1.0, 1.0, 1.0, 1.0);
        gl.enable(gl.DEPTH_TEST);
	    //gl.enable(gl.CULL_FACE);
		
		var $=this;
		$.vMatrix = mat4.create();
		$.mMatrix = mat4.create();
		$.pMatrix = mat4.create();
		$.mvMatrix = mat4.create();
		$.tempMatrix = mat4.create();
		$.normalMatrix = mat3.create();
		$.camView = null;
		$.camAngle = 0;

		this.targetPos = vec4.create();
	}
	
	
	/*used to resize the canvas*/
	p.resize = function(w,h)
	{
		var $=this;
		var c=$.canvas;
		var ctx=$.ctx;
	
		c.width = w;
		c.height = h;
		c.style.width = c.width+'px';
		c.style.height = c.height+'px';
		$.w=c.width;
		$.h=c.height;
		
		/*hacky way to reset projection matrix*/
		this.resetView(false);

	}	

	 p.initGL = function() 
	{
		try {
            this.gl = this.canvas.getContext("webgl");
			//this.gl = WebGLDebugUtils.makeDebugContext(this.canvas.getContext("webgl"));            
        } catch (e) {
        }
        if (!this.gl) {
            alert("Could not initialise WebGL");
        }
    }

    p.getShader = function(gl, id) 
	{
        var shaderScript = document.getElementById(id);
        if (!shaderScript) {
			console.error("Failed to find shader "+id);
            return null;
        }

        var str = "";
        var k = shaderScript.firstChild;
        while (k) {
            if (k.nodeType == 3) {
                str += k.textContent;
            }
            k = k.nextSibling;
        }

        var shader;
        if (shaderScript.type == "x-shader/x-fragment") {
            shader = gl.createShader(gl.FRAGMENT_SHADER);
        } else if (shaderScript.type == "x-shader/x-vertex") {
            shader = gl.createShader(gl.VERTEX_SHADER);
        } else {
            return null;
        }

        gl.shaderSource(shader, str);
        gl.compileShader(shader);

        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            alert(gl.getShaderInfoLog(shader));
            return null;
        }

        return shader;
    }

    var shaderProgram;

	p.initShaders = function() 
	{
		var $=this;
		var gl = $.gl;
		var fragmentShader = $.getShader(gl, "shader-fs");
		var vertexShader = $.getShader(gl, "shader-vs");

		var shaderProgram = gl.createProgram();
		gl.attachShader(shaderProgram, vertexShader);
		gl.attachShader(shaderProgram, fragmentShader);
		gl.linkProgram(shaderProgram);

		if (!gl.getProgramParameter(shaderProgram, gl.LINK_STATUS)) {
			alert("Could not initialise shaders");
		}

		gl.useProgram(shaderProgram);
			
		shaderProgram.vertexPositionAttribute = gl.getAttribLocation(shaderProgram, "aPosition");
		gl.enableVertexAttribArray(shaderProgram.vertexPositionAttribute);

		shaderProgram.vertexColorAttribute = gl.getAttribLocation(shaderProgram, "aColor");
		gl.enableVertexAttribArray(shaderProgram.vertexColorAttribute);
	   
		shaderProgram.vertexNormalAttribute = gl.getAttribLocation(shaderProgram, "aNormal");
		gl.enableVertexAttribArray(shaderProgram.vertexNormalAttribute);

		shaderProgram.uMatrix = gl.getUniformLocation(shaderProgram, "uMatrix");        

		/*save*/
		$.shaderProgram = shaderProgram;
	}
    
	/*resets viex matrixes*/
	p.resetView = function(resetCam)
	{
		/*compute bounding box over all components*/
		var components = results_3d.components;
		
		/*find first component with defined box*/
		var c=-1;
		for (var i=0;i<components.length;i++) if (components[i].box1) {c=i;break;}
		if (c<0) return;
		
		var box1 = components[c].box1;
		var box2 = components[c].box2;
		
		for (var i=c+1;i<components.length;i++)
		{
			var comp = components[i];
			for (var d=0;d<3;d++)
			{
				if (comp.box1[d]<box1[d]) box1[d]=comp.box1[d];
				if (comp.box2[d]>box2[d]) box2[d]=comp.box2[d];
			}
		}
		
		/*range*/
		var xr = box2[0]-box1[0];
		var yr = box2[1]-box1[1];
		var zr = box2[2]-box1[2];

		/*model view matrix*/
		mat4.identity(this.mMatrix);

		if (resetCam)
		{
			this.camAngle = 0*Math.PI/180.0;
			var r=25;
			var x0=[box1[0]+0.5*xr,box1[1]+0.5*yr,box1[2]+0.5*zr];

			var camMatrix = mat4.create();
			mat4.lookAt(camMatrix, [x0[0], x0[1], x0[2]-r], [x0[0],x0[1],x0[2]], [0,1,0]);
			mat4.invert(this.vMatrix,camMatrix);
		}
		
		//mat4.perspective(this.pMatrix,90*Math.PI/180.0, width / height, 0.1, 100);
		var ar = this.w/this.h;
		var x1,x2,y1,y2,z1,z2;
		xr*=10;
		yr*=10;
		zr*=10;
		x1 = ar*(box1[0]-xr); x2 = ar*(box2[0]+xr);
		y1 = box1[1]-yr; y2 = box2[1]+yr;
		z1 = box1[2]-zr; z2 = 3*zr;
		mat4.ortho(this.pMatrix,x1,x2,y1,y2,z1,z2);
		/*console.log(x1+" "+x2);
		console.log(y1+" "+y2);
		console.log(z1+" "+z2);
		*/
					   
		/*set light*/
		var gl=this.gl;
		gl.uniform3f(gl.getUniformLocation(this.shaderProgram,"uLightPosition"), -10,-10,0);
		gl.uniform3f(gl.getUniformLocation(this.shaderProgram,"uLightColor"), 1.0, 1.0, 1.0);
		gl.uniform3f(gl.getUniformLocation(this.shaderProgram,"uAmbientLight"), 1.0, 1.0, 1.0);
		
		this.drawScene();
	}
	
	/*drawScene ***************************************
	This function sets the model-view matrix and then calls drawComponent to actually
	plot the data. The idea is that if we have multiple objects, we could loop over them by
	calling drawComponent in a loop	*/
	p.drawScene = function() 
	{
		/*return if no data*/
		if (results_3d.components.length<1) return;
		
		var $=this;
		var gl = $.gl;
		var width = gl.drawingBufferWidth;
		var height = gl.drawingBufferHeight;
		gl.viewport(0, 0, width, height);	/*draw in entire context*/
		
		/*reset depth bit*/
		gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

		/*U=U*P, need to multiply from left due to column major?*/
	
		/*premultiply matrix as it is used by all vertices
		pMatrix * vMatrix * mMatrix*/	
		mat4.mul($.tempMatrix,$.vMatrix,$.mMatrix);
		mat4.mul($.mvMatrix,$.pMatrix,$.tempMatrix);
		
		var shaderProgram = $.shaderProgram;
		
		gl.uniformMatrix4fv(gl.getUniformLocation(shaderProgram,"mvMatrix"), false, $.mvMatrix);
		
		/*set normals*/
		mat3.normalFromMat4($.normalMatrix,$.mMatrix);
		gl.uniformMatrix3fv(gl.getUniformLocation(shaderProgram,"uNormalMatrix"), false, $.normalMatrix);
		
		/*this is what will actually do the drawing*/
		$.drawComponent(0);
		
		var x= [0, 0, 0, 0];
		var str=mat4.str($.mvMatrix);			
		//console.log(str);
	}

	/*sets active component by copying geometry to the GPU*/
	p.drawComponent = function(comp_index)
	{
		var comp = results_3d.components[comp_index];
		if (!comp.visible) return;
			
		var gl = this.gl;
		var shaderProgram = this.shaderProgram;

		gl.bindBuffer(gl.ARRAY_BUFFER, comp.positionBuffer);
		gl.vertexAttribPointer(shaderProgram.vertexPositionAttribute, 3, gl.FLOAT, false, 0, 0);
		
		/*normal vectors*/
		gl.bindBuffer(gl.ARRAY_BUFFER, comp.normalsBuffer);
		gl.vertexAttribPointer(shaderProgram.vertexNormalAttribute, 3, gl.FLOAT, false, 0, 0);
		
		/*draw elements, types include gl.TRIANGLES, POINTS, LINES*/
		if (comp.draw_mode&gl.LINES)
		{
			gl.uniform1i(gl.getUniformLocation(shaderProgram, "black"), 1);
			gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, comp.lineIndexBuffer);	
			gl.drawElements(gl.LINES, comp.n_lines*2, gl.UNSIGNED_SHORT,0);		
		}
		
		if (comp.draw_mode&gl.TRIANGLES)
		{
			gl.uniform1i(gl.getUniformLocation(shaderProgram, "black"), 0);
			gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, comp.triIndexBuffer);	
			gl.drawElements(gl.TRIANGLES, comp.n_elements*3, gl.UNSIGNED_SHORT,0);
			//gl.drawElements(gl.POINTS,comp.triIndexBuffer.numItems*3,gl.UNSIGNED_SHORT,0);
		}
	
	}
	
	
	/*copies geometry data to the GPU*/
	p.setComponent = function(data,comp) 
	{        
		/*append new component if not specified*/
		if (!comp)
		{
			comp = {};
			results_3d.components.push(comp);
		}
		
		var vertices = data.vertices;
		var indexes = data.indexes;
		var normals = data.normals;
		var scalars = data.scalars;
		
		var gl = this.gl;
		var shaderProgram = this.shaderProgram;
		
		comp.visible = true;
		if (data.draw_mode) comp.draw_mode = data.draw_mode;
		else comp.draw_mode = gl.TRIANGLES|gl.LINES;
		
		comp.vertices=vertices;
		comp.n_nodes = vertices.length/3;
		comp.n_elements = indexes.length/3;
		comp.positionBuffer = gl.createBuffer();
		gl.bindBuffer(gl.ARRAY_BUFFER, comp.positionBuffer);
		gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);
		
		/*create element index buffer*/
		comp.triIndexBuffer = gl.createBuffer();
		gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, comp.triIndexBuffer);
		gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(indexes), gl.STATIC_DRAW);
			
		/*index buffer for lines*/
		var line_indexes = [];
		for (var i=0;i<comp.n_elements;i++)
		{
			line_indexes.push(indexes[i*3],indexes[i*3+1],
							  indexes[i*3+1],indexes[i*3+2],
							  indexes[i*3+2],indexes[i*3]);
		}
		comp.lineIndexBuffer = gl.createBuffer();
		gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, comp.lineIndexBuffer);
		gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(line_indexes), gl.STATIC_DRAW);
		comp.n_lines = comp.n_elements*3;
		
		/*normal vectors*/
		comp.normalsBuffer = gl.createBuffer();
		gl.bindBuffer(gl.ARRAY_BUFFER, comp.normalsBuffer);
		gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(normals), gl.STATIC_DRAW);
		gl.vertexAttribPointer(shaderProgram.vertexNormalAttribute, 3, gl.FLOAT, false, 0, 0);
	
		/*set scalars*/
		if (scalars)
			comp.scalars = scalars;
		else 
			comp.scalars = [];
		
		/*storage for scalars, initialize to black*/
		comp.colorBuffer = gl.createBuffer();
		comp.colorBuffer.itemSize = 3;
		comp.colorBuffer.numItems = data.length;  
		gl.bindBuffer(gl.ARRAY_BUFFER, comp.colorBuffer);
		var colors = [];
		for (var i=0;i<comp.n_nodes;i++) {colors.push(0.5);colors.push(0.5);colors.push(0.5);}
		gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(colors), gl.STATIC_DRAW);
		gl.vertexAttribPointer(shaderProgram.vertexColorAttribute, comp.colorBuffer.itemSize, gl.FLOAT, false, 0, 0);

		
		/*get bounding box*/
		var box1=[vertices[0],vertices[1],vertices[2]];
		var box2=[box1[0],box1[1],box1[2]];
		for (var i=0;i<vertices.length;i+=3)
		{
			for (var j=0;j<3;j++)
			{
				if (box1[j]>vertices[i+j]) box1[j]=vertices[i+j];
				if (box2[j]<vertices[i+j]) box2[j]=vertices[i+j];
			}
		}
		comp.box1 = box1;
		comp.box2 = box2;
		
	}

	/*converts scalar data to RGB*/
	p.convertScalars = function(var_index)
	{
		/*only the first component contains scalar data*/
		var comp = results_3d.components[0];
		
		if (comp.scalars.length==0) return;

		var scalars = comp.scalars[var_index];
		var data = scalars.data;

		var colors = []
		var range=scalars.max-scalars.min;
		if (!range) range=1;
		
		for (var i=0; i < data.length; i++) 
		{
			var f=(data[i]-scalars.min)/range;
			if (f>1.0) f=1.0;
			else if (f<0.0) f=0.0;
			f = (Math.floor(f*scalars.levels))/scalars.levels;
			var a=(1-f)/0.25;	//invert and group
			var X=Math.floor(a);	//this is the integer part
			var Y=Math.floor(255*(a-X)); //fractional part from 0 to 255
			switch(X)
			{
				case 0: r=255;g=Y;b=0;break;
				case 1: r=255-Y;g=255;b=0;break;
				case 2: r=0;g=255;b=Y;break;
				case 3: r=0;g=255-Y;b=255;break;
				case 4: r=0;g=0;b=255;break;
				default: r=0;g=0;b=0;break;
			}
			
			colors.push(r/255.0);
			colors.push(g/255.0);
			colors.push(b/255.0);
		}
		
		/*colors*/
		var gl = this.gl;
		var shaderProgram = this.shaderProgram;
		gl.bindBuffer(gl.ARRAY_BUFFER, comp.colorBuffer);
		gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(colors), gl.STATIC_DRAW);
		gl.vertexAttribPointer(shaderProgram.vertexColorAttribute, comp.colorBuffer.itemSize, gl.FLOAT, false, 0, 0);	
	}
		
 
	p.handleMouseDown = function(event)
	{
		this.mouseDown = true;
		this.lastMouseX = null;
		this.lastMouseY = null;
		this.mouseClick = null;

		if (event.buttons==1) this.mouseClick = CLICK.LEFT;
		else if (event.buttons==2) this.mouseClick = CLICK.RIGHT;
		else if (event.buttons==3 || event.buttons==4) this.mouseClick = CLICK.MIDDLE;
		this.lastMouseX = event.clientX;
		this.lastMouseY = event.clientY;
	}

	p.handleMouseUp = function (event) 
	{
		this.mouseDown = false;
	}

    p.handleMouseMove = function(event) 
	{
		if (!this.mouseDown)  return;
		
		event.preventDefault();
	    var newX = parseInt(event.clientX);
		var newY = parseInt(event.clientY);

		var deltaX = newX - this.lastMouseX;
		var deltaY = newY - this.lastMouseY;
    
		if (this.mouseClick == CLICK.LEFT)
		{
			/*rotation with left click*/
			mat4.rotateY(this.vMatrix, this.vMatrix, degToRad(deltaX / 2.0));
			mat4.rotateX(this.vMatrix, this.vMatrix, degToRad(deltaY / 2.0));
		}
		else if (this.mouseClick == CLICK.RIGHT)
		{
			/*need to return to avoid setting lastMouseY*/
			if (Math.abs(deltaY)<2) return;
			
			f=1.1;
			if (deltaY<0) f=1.0/f;
			mat4.scale(this.vMatrix, this.vMatrix, [f,f,f]);
			
		}
		else if (this.mouseClick == CLICK.MIDDLE)
		{
			/*translation*/
			mat4.translate(this.vMatrix, this.vMatrix, [-deltaX/15.0,-deltaY/15.0,0]);
		}

		
		this.lastMouseX = newX
		this.lastMouseY = newY;
		
	}

  
	this.init(id);

}	/*CANVAS*/

/**********************************************************************************/



/*4x4 matrix vec4 mult
0 4 8 12
1 5 9 13
2 6 10 14
3 7 11 15
*/
/*matrix math function*/
mat4.multiplyVec4= function (r,M,x)
{
	if (!x) x=r;
	var t=[0,0,0,0];
	for (var i=0;i<4;i++)
	{
		t[i] = M[i]*x[0] + M[i+4]*x[1] + 
		       M[i+8]*x[2] + M[i+12]*x[3];
	}
	
	/*copy*/
	for (var i=0;i<4;i++) r[i]=t[i];
}

