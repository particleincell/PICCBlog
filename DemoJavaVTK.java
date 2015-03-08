import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JToggleButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import vtk.vtkNativeLibrary;
import vtk.vtkPanel;
import vtk.vtkActor;
import vtk.vtkSphere;
import vtk.vtkSphereSource;
import vtk.vtkSampleFunction;
import vtk.vtkContourFilter;
import vtk.vtkPlane;
import vtk.vtkCutter;
import vtk.vtkLookupTable;
import vtk.vtkPolyDataMapper;

/* ************************************************************
 * Demo applications showcasing how to use VTK with Java
 * 
 * Based on SimpleVTK.java example distributed with VTK
 * 
 * For more information see:
 * http://www.particleincell.com/2011/vtk-java-visualization
 * 
 * Information about VTK can be found at:
 * http://vtk.org/
 * 
 * ***********************************************************/

public class DemoJavaVTK extends JPanel implements ActionListener 
{
    private static final long serialVersionUID = 1L;
    private vtkPanel renWin;
    private vtkActor cutActor;
    private vtkActor isoActor;
	
    private JPanel buttons;
    private JToggleButton slicesButton;
    private JToggleButton isoButton;
    private JButton exitButton;

    /* Load VTK shared librarires (.dll) on startup, print message if not found */
    static 
    {				
        if (!vtkNativeLibrary.LoadAllNativeLibraries()) 
	{
	       for (vtkNativeLibrary lib : vtkNativeLibrary.values()) 
		{
                	if (!lib.IsLoaded()) 
				System.out.println(lib.GetLibraryName() + " not loaded");    
		}
			
		System.out.println("Make sure the search path is correct: ");
		System.out.println(System.getProperty("java.library.path"));
        }
        vtkNativeLibrary.DisableOutputWindow(null);
    }

    /* Constructor - generates visualization pipeline and adds actors*/
    public DemoJavaVTK() 
    {
        super(new BorderLayout()); /* large center and small border areas*/
		
        double radius = 0.8;		/*sphere radius*/
		
	/**** 1) INPUT DATA: Sphere Implicit Function ****/
	vtkSphere sphere = new vtkSphere();
	sphere.SetRadius(radius);
		
	vtkSampleFunction sample = new vtkSampleFunction();
	sample.SetSampleDimensions(50,50,50);
	sample.SetImplicitFunction(sphere);
		
	/**** 2) PIPELINE 1: Isosurface Actor ****/
		
	/* contour filter - will generate isosurfaces from 3D data*/
	vtkContourFilter contour = new vtkContourFilter();
	contour.SetInputConnection(sample.GetOutputPort());
	contour.GenerateValues(3,0,1);
		
	/* mapper, translates polygonal representation to graphics primitives */
	vtkPolyDataMapper isoMapper = new vtkPolyDataMapper();
        isoMapper.SetInputConnection(contour.GetOutputPort());
		
	/*isosurface actor*/
        isoActor = new vtkActor();
        isoActor.SetMapper(isoMapper);
	
	/**** 3) PIPELINE 2: Cutting Plane Actor ****/
		
	/* define a plane in x-y plane and passing through the origin*/
	vtkPlane plane = new vtkPlane();
	plane.SetOrigin(0,0,0);
	plane.SetNormal(0,0,1);
		
	/* cutter, basically interpolates source data onto the plane */
	vtkCutter planeCut = new vtkCutter();
	planeCut.SetInputConnection(sample.GetOutputPort());
	planeCut.SetCutFunction(plane);
	/*this will actually create 3 planes at the subspace where the implicit
	 * function evaluates to -0.7, 0, 0.7 (0 would be original plane). In 
	 * our case this will create three x-y planes passing through 
	 * z=-0.7, z=0, and z=+0.7*/
	planeCut.GenerateValues(3,-0.7,0.7);
		
	/* look up table, we want to reduce number of values to get discrete bands */
	vtkLookupTable lut = new vtkLookupTable();
	lut.SetNumberOfTableValues(5);
		
	/* mapper, using our custom LUT */
	vtkPolyDataMapper cutMapper = new vtkPolyDataMapper();
        cutMapper.SetInputConnection(planeCut.GetOutputPort());
	cutMapper.SetLookupTable(lut);
		
	/* cutting plane actor, looks much better with flat shading */
	cutActor = new vtkActor();
        cutActor.SetMapper(cutMapper);
	cutActor.GetProperty().SetInterpolationToFlat();
		
	/**** 4) PIPELINE 3: Surface Geometry Actor ****/
		
	/* create polygonal representation of a sphere */
	vtkSphereSource surf = new vtkSphereSource();
	surf.SetRadius(radius);
		
	/* another mapper*/
	vtkPolyDataMapper surfMapper = new vtkPolyDataMapper();
	surfMapper.SetInputConnection(surf.GetOutputPort());
		
	/* surface geometry actor, turn on edges and apply flat shading*/
	vtkActor surfActor = new vtkActor();
	surfActor.SetMapper(surfMapper);
	surfActor.GetProperty().EdgeVisibilityOn();
	surfActor.GetProperty().SetEdgeColor(0.2,0.2,0.2);
	surfActor.GetProperty().SetInterpolationToFlat();

	/**** 5) RENDER WINDOW ****/
		
	/* vtkPanel - this is the interface between Java and VTK */
	renWin = new vtkPanel();
		
	/* add the surface geometry plus the isosurface */
	renWin.GetRenderer().AddActor(surfActor);
	renWin.GetRenderer().AddActor(isoActor);
		
	/* the default zoom is whacky, zoom out to see the whole domain */
        renWin.GetRenderer().GetActiveCamera().Dolly(0.2); 
	renWin.GetRenderer().SetBackground(1, 1, 1);
		
	/**** 6) CREATE PANEL FOR BUTTONS ****/
	buttons  = new JPanel();
	buttons.setLayout(new GridLayout(1,0));
		
        /* isosurface button, clicked by default */
	isoButton = new JToggleButton("Isosurfaces",true);
        isoButton.addActionListener(this);
		
	/* cutting planes button */
        slicesButton = new JToggleButton("Slices");
        slicesButton.addActionListener(this);
		
	/* exit button */
	exitButton = new JButton("Exit");
        exitButton.addActionListener(this);
		
	/* add buttons to the panel */
	buttons.add(isoButton); 
	buttons.add(slicesButton);
	buttons.add(exitButton); 

	/**** 7) POPULATE MAIN PANEL ****/
        add(renWin, BorderLayout.CENTER);
        add(buttons, BorderLayout.SOUTH);	
    }

    /* ActionListener that responds to button clicks
     * Toggling iso/slices buttons results in addition or removal
     * of the corresponding actor */
    public void actionPerformed(ActionEvent e) 
    {
	/*cutting planes button, add or remove cutActor */
	if (e.getSource().equals(slicesButton))
	{
		if (slicesButton.isSelected())
			renWin.GetRenderer().AddActor(cutActor);
		else
			renWin.GetRenderer().RemoveActor(cutActor);
			
		renWin.Render();
	}
	/*isosurface button, add or remove isoActor */
	else if (e.getSource().equals(isoButton))
	{
		if (isoButton.isSelected())
			renWin.GetRenderer().AddActor(isoActor);
		else
			renWin.GetRenderer().RemoveActor(isoActor);
		renWin.Render();
	}
	/*exit button, end application */
	else if (e.getSource().equals(exitButton)) 
	{
            System.exit(0);
        }
    }

    /* main, creates a new JFrame and populates it with the DemoJavaVTK panel */
    public static void main(String s[]) 
    {
        SwingUtilities.invokeLater(new Runnable() 
	{
            @Override
            public void run() 
	    {
                JFrame frame = new JFrame("Java and VTK Demo");
                frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                frame.getContentPane().setLayout(new BorderLayout());
                frame.getContentPane().add(new DemoJavaVTK(), BorderLayout.CENTER);
                frame.setSize(400, 400);
                frame.setLocationRelativeTo(null);
                frame.setVisible(true);
            }
        });
    }
}
