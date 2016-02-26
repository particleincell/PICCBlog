/*support for file reading*/
  function File() {
	var t=this;
	var p=File.prototype;  
  }
 
/*reads in .xml file*/
function loadXML(filename)
{
	if (!filename) filename=document.getElementById("file_name").value;
	
	localStorage.file_name = filename;
	
	var xhttp;
	
	if (window.XMLHttpRequest)
	{
	  xhttp=new XMLHttpRequest();
	}
	else {alert("window.XMLHttpRequest not supported!");return;}
	 
	xhttp.onload = function (e) {
	if (xhttp.readyState === 4) {
    if (xhttp.status === 200) {
		var xml = xhttp.responseXML;
		parseVTK(xml);				
    } else {
      log(xhttp.statusText);
    }
  }
};
	xhttp.onerror = function (e) {
	log(xhttp.statusText);};
	xhttp.open("GET", filename, true);
	xhttp.send(null);
	
}

/*******VTK LOADER********************/
function parseVTK(xml)
{
	if (typeof(xml.firstChild.children)=="undefined") {log_error("Input file does not seem to be in XML format");return false;}
		
	var grid = xml.firstChild.children[0];
	var file_type;
	if (grid.tagName=="UnstructuredGrid") file_type="UG";
	else if (grid.tagName=="PolyData") file_type="PD";
	else {log_error("Expected &lt;UnstructuredGrid&gt; or &lt;PolyData&gt; instead found &lt;"+grid.tagName+"&gt;");return false;}
	
	var vertices = [];
    var node_indexes=[];
    var normals = [];
    var indexes = [];
	var scalars = [];
	    
	var pieces = xml.getElementsByTagName("Piece");
	
	for (var p=0;p<pieces.length;p++)
	{
		//if (p!=2) continue;
		var piece = pieces[p];
		var n_nodes = piece.getAttribute("NumberOfPoints");
		var n_cells;
		
		var cell0 = indexes.length/3;	/*number of cells prior to this pieces*/
		var node0 = vertices.length/3; /*number of nodes prior to this piece*/
		
		if (file_type=="UG") n_cells=piece.getAttribute("NumberOfCells");
		else n_cells=piece.getAttribute("NumberOfPolys");
		
		console.log("Loading piece "+p+" with "+n_nodes+" nodes and "+n_cells+" cells");
		/*read points*/
		var points = piece.getElementsByTagName("Points")[0];
		var da = points.getElementsByTagName("DataArray")[0];
	 
		var buff = (da.innerHTML).trim();
		
		/*split by spaces*/
		data = buff.split(/\s+/);
        
		for (var n=0;n<n_nodes;n++)
		{
            /*read node data*/
			var x = parseFloat(data[n*3]);
            var y = parseFloat(data[n*3+1]);
            var z = parseFloat(data[n*3+2]);
                        
            vertices.push(x);
            vertices.push(y);
            vertices.push(z);
		}
		
		/*unstructured grid*/
		if (file_type=="UG")
		{
			/*read elements*/
			var cells = piece.getElementsByTagName("Cells")[0];
			var da = cells.getElementsByTagName("DataArray")[0];
			var buff = (da.innerHTML).trim();
			/*split by spaces*/
			data = buff.split(/\s+/);
			
			for (var c=0;c<n_cells*3;c++)
			{
				var n = parseInt(data[c]);
				indexes.push(n+node0);	
			}
		}
		else
		{
			/*polydata*/
			/*read elements*/
			var polys = piece.getElementsByTagName("Polys")[0];
			var da = polys.getElementsByTagName("DataArray");
			
			var connect;
			var offsets;
			for (var i=0;i<da.length;i++)
			{
				var name = da[i].getAttribute("Name");
				var buff = (da[i]).innerHTML.trim();
				buff = buff.split(/\s+/);
				if (name=="connectivity") connect = buff;
				else if (name=="offsets") offsets = buff;
			}
			
			var last = 0;
			for (var i=0;i<n_cells;i++)
			{
				var cur = parseInt(offsets[i]);
				if ((cur-last)!=3) {log_error("Only triangles are supported, ("+(cur-last)+")");continue;}
				
				for (var j=0;j<3;j++)
				{
					var n = parseInt(connect[last+j]);
					indexes.push(n+node0);					
				}
				last = cur;
			}
		}
		
		/*save user view data*/
		
		
		/*load point data*/
		var point_data = piece.getElementsByTagName("PointData")[0];
		if (point_data)
		for (var i=0;i<point_data.children.length;i++)
		{
			var da = point_data.children[i];
			var nc = parseInt(da.getAttribute("NumberOfComponents"));
			if (nc>1) continue;
			var name = da.getAttribute("Name");
			var buff = (da.innerHTML).trim();
			var sdata = buff.split(/\s+/);
			/*do we already have this scalar?*/
			
			var data=[];
			for (var j=0;j<scalars.length;j++) {if (scalars[j].name==name) {scalars[j].data=[];data=scalars[j].data;break;}}
			/*if didn't find it*/
			if (data.length==0) {scalars.push({name:name,data:[]}); data = scalars[scalars.length-1].data;}
			
			for (var n=0;n<n_nodes;n++)
			{
				data.push(parseFloat(sdata[n]));		
			}					
		}
					
		/*now load cell data, converted to point data*/
		var cell_data = piece.getElementsByTagName("CellData")[0];
		var data;
		for (var i=0;i<cell_data.children.length;i++)
		{
			var da = cell_data.children[i];
			var nc = parseInt(da.getAttribute("NumberOfComponents"));
			if (nc>1) continue;
			var name = da.getAttribute("Name")+"*";
			var buff = (da.innerHTML).trim();
			var sdata = buff.split(/\s+/);
			var c_scalars=[];
		
			for (var c=0;c<n_cells;c++) 
				c_scalars.push(parseFloat(sdata[c]));						
			
			/*our buffer needs point data, convert*/
			var counts = [];
			
			var data = [];
			for (var j=0;j<scalars.length;j++) {if (scalars[j].name==name) {data=scalars[j].data;break;}}
			/*if didn't find it*/
			if (data.length==0) {scalars.push({name:name,data:[]}); data = scalars[scalars.length-1].data;}
			
			/*add n_nodes zeros to the end of data*/
			for (var n=0;n<n_nodes;n++) {data.push(0);counts.push(0);}
			for (var c=0;c<n_cells;c++) {
				for (var j=0;j<3;j++)
				{
					var n = indexes[(cell0+c)*3+j];
					data[n]+= c_scalars[c];
					counts[n-node0]++;
				}
			}
			for (var n=0;n<n_nodes;n++) {if (counts[n]>0) data[node0+n]/=counts[n]; else data[node0+n]=0;}					
		}
		
	}	/*end of a piece*/
	
	/*make sure we got scalars*/
	if (scalars.length==0)
	{
		console.error("No scalars found, setting to zero");
		scalars[0].name="empty";
		for (var i=0;i<vertices.length;i++) scalars[0].data.push(0);		
	}
	
	/*compute normals*/
	var c_normals = [];
	n_cells = indexes.length/3;
	n_nodes = vertices.length/3;
	for (var c=0;c<n_cells;c++)
	{
		var n1 = indexes[c*3];
		var n2 = indexes[c*3+1];
		var n3 = indexes[c*3+2];

		var x1 = [vertices[n1*3],vertices[n1*3+1],vertices[n1*3+2]];
		var x2 = [vertices[n2*3],vertices[n2*3+1],vertices[n2*3+2]];
		var x3 = [vertices[n3*3],vertices[n3*3+1],vertices[n3*3+2]];
		
		var a = [0,0,0];
		var b = [0,0,0];
		for (var dim=0;dim<3;dim++) {a[dim]=x2[dim]-x1[dim];b[dim]=x3[dim]-x1[dim];}
		
		c_normals.push(a[1]*b[2]-a[2]*b[1]);
		c_normals.push(-(a[0]*b[2]-a[2]*b[0]));
		c_normals.push(a[0]*b[1]-a[1]*b[0]);				
	}
	
	//for (var c=0;c<n_cells;c++) {c_normals[c*3]=1;c_normals[c*3+1]=0;c_normals[c*3+2]=0;}
	/*convert to node data*/
	var counts = [];
	for (var n=0;n<n_nodes;n++) {normals.push(0);normals.push(0);normals.push(0);counts.push(0);}
	for (var c=0;c<n_cells;c++) 
	{
		var n1 = indexes[c*3];
		var n2 = indexes[c*3+1];
		var n3 = indexes[c*3+2];
		for (var dim=0;dim<3;dim++)	/*add xyz*/
		{
			normals[n1*3+dim]+= c_normals[c*3+dim];
			normals[n2*3+dim]+= c_normals[c*3+dim];
			normals[n3*3+dim]+= c_normals[c*3+dim];
		}
		counts[n1]++;
		counts[n2]++;
		counts[n3]++;		
	}
	
	for (var n=0;n<n_nodes;n++) {if (counts[n]>0) for (var dim=0;dim<3;dim++) normals[n*3+dim]/=counts[n];}
	
	/*normalize*/
	for (var n=0;n<n_nodes;n++) {
		var j=n*3; 
		var nx=normals[j];
		var ny=normals[j+1];
		var nz=normals[j+2];
		var mag=Math.sqrt(nx*nx+ny*ny+nz*nz);
		
		for (var dim=0;dim<3;dim++) normals[j+dim]/=mag;
	}
	
	/*set scalar ranges*/
	for (var i=0;i<scalars.length;i++)
	{
		var s = scalars[i];
		/*get data range*/
		s.range_min = s.data[0];
		s.range_max = s.data[0];
		for (var j=0;j<data.length;j++)
		{
			if (s.data[j]<s.range_min) s.range_min=s.data[j];
			if (s.data[j]>s.range_max) s.range_max=s.data[j];
		}
		s.min = s.range_min;
		s.max = s.range_max;
		s.levels = 32;
		s.log = false;
	}
	/*sanity check*/
	for (var i=0;i<indexes.lengths;i++)
		if (indexes[i]<0 || indexes[i]>=vertices.length) console.log("Bad Index: "+indexes[i]);
	
	data = {vertices,node_indexes,normals,indexes,scalars};
	
	/*update results with this data*/
	results_3d.setData(data);

	return true;
}



/*------- FILE READING ------*/
  function handleFileSelect(evt) {
    evt.stopPropagation();
    evt.preventDefault();
	log_clear();
	var files = evt.dataTransfer.files; // FileList object.
	evt.target.classList.remove("ondrag");
	
    /*read only first file*/
	f=files[0];
	
	/*parse file*/
	var reader = new FileReader();
	
	// Closure to capture the file information.
	reader.onload = (function(file,target) {
		return function(e) {
		if (target.id=="drop")
		{
			/*this is not going to work, need to convert target.result from text to XML*/
			var parser = new DOMParser();
			var xml = parser.parseFromString(e.target.result, "application/xml")
			if (parseVTK(xml))
			log("Done reading "+file.name);		        
		}		
		};
      })(f,evt.target);
	  
	var data = reader.readAsText(f);    
  }

    
  function handleDragOver(evt) {
    evt.dataTransfer.dropEffect = 'copy'; 
	evt.stopPropagation();
    evt.preventDefault();
	evt.target.classList.add("ondrag");
	}

	function handleDragLeave(evt) {
	/*evt.target seems to be the new element over the mouse, currentTarget is the original*/
	evt.currentTarget.classList.remove("ondrag");
		
	}
