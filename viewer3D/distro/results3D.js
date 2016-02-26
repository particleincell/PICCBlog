
var results_3d;

/**** INIT ***************/
function init3D(tab_index)
{
	var div = insertTab(tab_index); 
	results_3d = new results_3d(div, tab_index);
	tab_data[tab_index] = results_3d;
}

 var CLICK = {
	  LEFT : {value: 0, name: "Left"}, 
	  RIGHT: {value: 1, name: "Right"}, 
	  MIDDLE : {value: 2, name: "Middle"}
	};

function degToRad(d) {return d*Math.PI/180.0;}
	
var results_3d = function(div, tab_index)
{
	var p = results_3d.prototype;
	
	p.init = function(div, tab_index)
	{
		/*copy content*/
		var source = document.getElementById("content3D");
		div.innerHTML=source.innerHTML;
		
		/*get file name from storage*/
		var file_name = getFromStorage("input_3d","surface.vtu");
			
		document.getElementById("c_3d").addEventListener("contextmenu", function(evt) {evt.preventDefault();});
		this.canvas = new Canvas3D("c_3d");
		
		document.getElementById("input_3d").value=file_name;
		document.getElementById("auto_range_3d").addEventListener("input", function(t){return function(evt) {t.changeAutoRange(evt.target);}}(this));
		
		document.getElementById("auto_reload_3d").addEventListener("change", function(t){return function(e){t.toggleAutoReload(e);}}(this));
		document.getElementById("reload_button_3d").addEventListener("click", reloadFile);    
		
		/*hooks for variable data*/
		document.getElementById("vars_3d").addEventListener("change", function(t){return function(evt){var index=evt.target.selectedIndex;t.changeVar(index);}}(this));
		var f_range = function(t){return function(e){t.changeRange(e);}}(this);
		document.getElementById("range_min_3d").addEventListener("input", f_range);
		document.getElementById("range_max_3d").addEventListener("input", f_range);
		document.getElementById("range_levels_3d").addEventListener("input", f_range);
		document.getElementById("range_log_3d").addEventListener("change", f_range);
		document.getElementById("domains_3d").addEventListener("change", function(t){return function(e){t.onDomainChange(e);}}(this));
		
		this.components = [{}];	/*first component is for geometry*/
		this.var_index = 0;	/*current variable*/
	}
	
	/*sets data to stuff from file*/
	p.setData = function(data)
	{
		this.canvas.setComponent(data,this.components[0]);
		
		this.initVars();
		this.canvas.convertScalars(this.var_index);
		this.canvas.resetView();
		
	}
	p.changeAutoRange = function(e)
	{
		var options = [5, 10, 30, 60, 2*60, 5*60, 10*60,15*60, 30*60,60*60];
		var val = parseInt(e.value);
		this.reload_interval = options[val-1];
		var hint;
		if (this.reload_interval<60) hint=this.reload_interval+" sec";
		else hint = (this.reload_interval/60)+" min";
		document.getElementById("rlegend_3d").innerHTML=hint;

		/*is auto on?*/
		if (document.getElementById("auto_reload_3d").checked) 
		{
			window.clearInterval(this.reload_interval_id);
			this.reload_interval_id = window.setInterval(reloadFile,this.reload_interval*1000);
		}
	}
 
	p.toggleAutoReload = function(evt)
	 {
		var e = evt.target;
		
		if (e.checked) this.reload_interval_id = window.setInterval(function(){
		reloadFile({target:document.getElementById("reload_button")})},
			this.reload_interval*1000);
		else window.clearInterval(this.reload_interval_id);
	 
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
  
	/*called on domain select change*/
	p.onDomainChange = function(evt)
	{
		for (var i=0;i<evt.target.options.length-1;i++)
		{
			var option = evt.target.options[i];
			this.components[i+1].visible=option.selected;
		}	
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

	/*called on tab activation*/
	p.refresh = function()
	{
		var l1 = this.components.length;
		
		/*generate boxes for domains*/
		var comps = [this.components[0]];	/*save geometry*/
		this.components = comps;
		
		var sel = document.getElementById("domains_3d");
		sel.innerHTML="";
			
		for (var i=0;i<world.domains.length;i++)
		{
			var dom = world.domains[i];
			var data = this.makeBox(dom.bound_lo.values, dom.bound_hi.values)
			this.canvas.setComponent(data);
			
			/*populate select box*/
			var option = document.createElement("option");
			option.text = dom.name.values;
			option.value = i+1;
			option.selected = true;
			sel.add(option);
			
		}
		
		var option = document.createElement("option");
		option.text = "";
		option.value = -1;
		//option.setAttribute("disabled","true");
		//option.setAttribute("hidden","true");
		sel.add(option)
		
		
		/*reset view if new domains added*/
		if (l1!=this.components.length) this.canvas.resetView();
		
		this.tick();		
	}
	
	
	/*generates points for a box*/
	p.makeBox = function(min, max)
	{
		var data = {};
		var x1 = parseFloat(min[0]);
		var y1 = parseFloat(min[1]);
		var z1 = parseFloat(min[2]);
		var x2 = parseFloat(max[0]);
		var y2 = parseFloat(max[1]);
		var z2 = parseFloat(max[2]);
		
		data.vertices = [x1,y1,z1,
						 x2,y1,z1,
						 x2,y2,z1,
						 x1,y2,z1,
						 x1,y1,z2,
						 x2,y1,z2,
						 x2,y2,z2,
						 x1,y2,z2];
		
		data.indexes = [1,4,0,
						1,4,5,
						5,1,2,
						6,5,2,
						2,3,6,
						3,7,6,
						0,4,3,
						3,4,7,
						0,3,1,
						1,3,2					
						];
		
		data.normals = [0,-1,0,
						0,-1,0,
						1,0,0,
						1,0,0,
						0,1,0,
						0,1,0,
						-1,0,0,
						-1,0,0,
						0,0,-1,
						0,0,-1,
						0,0,1,
						0,0,1];
		
		data.scalars = [[0,0,0,0,
						0,0,0,0]]

		data.draw_mode = this.canvas.gl.LINES;
		return data;
	}
	
	p.tick = function() 
	{
		if (activeTab==5)
		{
			requestAnimationFrame(function(t){return function(){t.tick();}}(this));
			this.checkSize();
			this.canvas.drawScene();
		}
	}

	this.init(div, tab_index);	
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
        
		var f_up=function(t){return function(e){if (activeTab==5) t.handleMouseUp(e);}}(this);
		var f_move=function(t){return function(e){if (activeTab==5) t.handleMouseMove(e);}}(this);
		
		this.canvas.addEventListener("mousedown",function(t){return function(e){t.handleMouseDown(e);}}(this));
		this.canvas.addEventListener("mouseup",f_up);
		this.canvas.addEventListener("mousemove",f_move);
		
		window.addEventListener("mouseup",f_up);
        window.addEventListener("mousemove", f_move);
		
		var gl = this.gl;
		gl.clearColor(1.0, 1.0, 1.0, 1.0);
        gl.enable(gl.DEPTH_TEST);
	    //gl.enable(gl.CULL_FACE);
		
		this.vMatrix = mat4.create();
		this.mMatrix = mat4.create();
		this.pMatrix = mat4.create();
		this.eyePos = vec4.create();
		this.targetPos = vec4.create();
	}
	
	p.clear = function() {
		this.ctx.fillStyle="#FFF";
		this.ctx.fillRect(0,0,this.w,this.h);
	}
	
	/*drawScene ***************************************/
	p.drawScene = function() 
	{
		/*return if no data*/
		if (results_3d.components.length<1) return;
		
		var gl = this.gl;
		var width = gl.drawingBufferWidth;
		var height = gl.drawingBufferHeight;
		gl.viewport(0, 0, width, height);	/*draw in entire context*/
		
		/*reset depth bit*/
		gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

		var camMatrix = mat4.create();
		mat4.lookAt(camMatrix, [this.eyePos[0], this.eyePos[1], this.eyePos[2]], 
					[this.targetPos[0],this.targetPos[1],this.targetPos[2]], this.lookUp);
		mat4.invert(this.vMatrix,camMatrix);

		/*perspective matrix*/
		mat4.identity(this.pMatrix);
		mat4.perspective(this.pMatrix,90*Math.PI/180.0, width / height, 0.1, 100);
	
		var shaderProgram = this.shaderProgram;
		gl.uniformMatrix4fv(gl.getUniformLocation(shaderProgram,"vMatrix"), false, this.vMatrix);
		gl.uniformMatrix4fv(gl.getUniformLocation(shaderProgram,"mMatrix"), false, this.mMatrix);
		gl.uniformMatrix4fv(gl.getUniformLocation(shaderProgram,"pMatrix"), false, this.pMatrix);

		/*reset black painting*/
		gl.uniform1i(gl.getUniformLocation(shaderProgram, "black"), 0);
		
		/*draw components*/
		for (var i=0;i<results_3d.components.length;i++)
		{
			this.drawComponent(i);			
		}
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
		
	//	$.ctx.translate(0,h);
	//	$.ctx.scale(1,-1);
	//	$.draw();	
	}	

	 p.initGL = function() 
	{
		try {
            //gl = canvas.getContext("webgl");
			this.gl = WebGLDebugUtils.makeDebugContext(this.canvas.getContext("webgl"));            
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
		shaderProgram.mMatrix = gl.getUniformLocation(shaderProgram, "mMatrix");
		
		/*save*/
		$.shaderProgram = shaderProgram;
	}
    
	/*resets viex matrixes*/
	p.resetView = function()
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
		
		console.log("Reset view bounding box");
		
		var xr = box2[0]-box1[0];
		var yr = box2[1]-box1[1];
		var zr = box2[2]-box1[2];
		
		/*range*/
		if (xr<0.01) xr=0.01;
		if (yr<0.01) yr=0.01;
		if (zr<0.01) zr=0.01;
	
		/*model view matrix*/
		mat4.identity(this.mMatrix);

		camAngle = 0*Math.PI/180.0;
		var r=1.5*Math.sqrt(xr*xr+yr*yr+zr*zr);
		var xc=[box1[0]+0.5*xr,box1[1]+0.5*yr,box1[2]+0.5*zr];

		console.log(xc);
		
		this.eyePos = vec4.fromValues(xc[0], xc[1], xc[2] - r * Math.cos(camAngle), 1);
		this.lookUp = vec4.fromValues(0,1,0,1);
		this.targetPos = vec4.fromValues(xc[0], xc[1], xc[2], 1);
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
		
		/*index buffer*/
		gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, comp.indexBuffer);
		
		/*normal vectors*/
		gl.bindBuffer(gl.ARRAY_BUFFER, comp.normalsBuffer);
		gl.vertexAttribPointer(shaderProgram.vertexNormalAttribute, 3, gl.FLOAT, false, 0, 0);
		
		var gl = this.gl;
		var shaderProgram = this.shaderProgram;

		/*draw elements, types include gl.TRIANGLES, POINTS, LINES*/
		if (comp.draw_mode&gl.LINES)
		{
			gl.uniform1i(gl.getUniformLocation(shaderProgram, "black"), 1);
			gl.drawElements(gl.LINES, comp.n_elements*3, gl.UNSIGNED_SHORT,0);		
		}
		
		if (comp.draw_mode&gl.TRIANGLES)
		{
			gl.uniform1i(gl.getUniformLocation(shaderProgram, "black"), 0);
			gl.drawElements(gl.TRIANGLES, comp.n_elements*3, gl.UNSIGNED_SHORT,0);		
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
		else comp.draw_mode = gl.TRIANGLES;
		
		comp.vertices=vertices;
		comp.n_nodes = vertices.length/3;
		comp.n_elements = indexes.length/3;
		comp.positionBuffer = gl.createBuffer();
		gl.bindBuffer(gl.ARRAY_BUFFER, comp.positionBuffer);
		gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);
		
		comp.indexBuffer = gl.createBuffer();
		gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, comp.indexBuffer);
		gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(indexes), gl.STATIC_DRAW);
		
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
		
	var mouseDown = false;
	var lastMouseX = null;
	var lastMouseY = null;
	var mouseClick = null;
	  
	p.handleMouseDown = function(event)
	{
		event.preventDefault();
		mouseDown = true;
		if (event.buttons==1) mouseClick = CLICK.LEFT;
		else if (event.buttons==2) mouseClick = CLICK.RIGHT;
		else if (event.buttons==3 || event.buttons==4) mouseClick = CLICK.MIDDLE;
		lastMouseX = event.clientX;
		lastMouseY = event.clientY;
	}

	p.handleMouseUp = function (event) 
	{
		mouseDown = false;
	}

    p.handleMouseMove = function(event) 
	{
		if (!mouseDown)  return;
		
		event.preventDefault();
	
		var newX = parseInt(event.clientX);
		var newY = parseInt(event.clientY);

		var deltaX = newX - lastMouseX;
		var deltaY = newY - lastMouseY;
		
		if (mouseClick == CLICK.LEFT)
		{
			var iMatrix = mat4.create();
			mat4.identity(iMatrix);
			/*rotation with left click*/
			mat4.rotateX(iMatrix, iMatrix, -degToRad(deltaY / 2.0));
			mat4.rotateY(iMatrix, iMatrix, -degToRad(deltaX / 2.0));
			mat4.multiplyVec4(this.eyePos, iMatrix, this.eyePos);
			
			/*repeat for lookup*/
			mat4.identity(iMatrix);
			mat4.rotateY(iMatrix, iMatrix, -degToRad(deltaX / 2.0));
			mat4.rotateX(iMatrix, iMatrix, -degToRad(deltaY / 2.0));
			mat4.multiplyVec4(this.lookUp, iMatrix, this.lookUp);
			
		}
		else if (mouseClick == CLICK.RIGHT)
		{
			/*need to return to avoid setting lastMouseY*/
			if (Math.abs(deltaY)<2) return;
			
			f=1.1;
			if (deltaY<0) f=1.0/f;
			var iMatrix = mat4.create();
			mat4.identity(iMatrix);
			mat4.scale(iMatrix, iMatrix, [f,f,f]);
			mat4.multiplyVec4(this.eyePos, iMatrix, this.eyePos);
			
		}
		else if (mouseClick == CLICK.MIDDLE)
		{
			/*translation*/
			var iMatrix = mat4.create();
			mat4.identity(iMatrix);
			mat4.translate(iMatrix, iMatrix, [deltaX/15.0,-deltaY/15.0,0]);
			mat4.multiplyVec4(this.eyePos, iMatrix, this.eyePos);
		}

	  
		
		lastMouseX = newX
		lastMouseY = newY;
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

/*parses the primary file*/
function load3D(buffer)
{
	data = parseVTK(buffer);
	results_3d.setData(data);

	
	/*switch for different file types here*/
	//results_xy.parseTable(buffer);		
}
