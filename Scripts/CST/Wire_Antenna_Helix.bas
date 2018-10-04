'#Language "WWB-COM"

Option Explicit

Sub Main

	SetUnit
	' Set up Frequency Sweep
	Dim Fstart As Double, Fend As Double, Fstep As Double
	Dim Fnum As Integer, cc As Integer, F() As String
	Fend = 130
	Fstart = 120
	Fstep = 1
	Fnum = Int(((Fend - Fstart)/Fstep) + 1)
	ReDim F(Fnum)
	SetSolver(F,Fstart,Fend,Fstep,Fnum,cc)
	' Set up Boundary Box
	Dim xminb, yminb, zminb, xmaxb, ymaxb, zmaxb
	xminb = 10
	yminb = 10
	zminb = 10
	xmaxb = 10
	ymaxb = 10
	zmaxb = 10
	' SetBoundBox (xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)
	' Create helix wire
	Dim h As String,ra As String,ri As String,phi As String,N As String,O As String, cPathStr As String, wPathStr As String
	h  = "10"
	ra = "5"
	ri = "5"
	phi = "30"
	N = "10"
	O = "1"
	cPathStr = "currentWire"
	wPathStr = "Wire"
	SetCoilWireCurrent(h,ra,ri,phi,N,O,cPathStr)
	' SetCoilWireCurrent(h,ra,ri,phi,N,O,wirePath)
	Dim circleStr As String, I As String
	' Set current along helix wire
	circleStr = "c"
	I = "1"
	SetCurrent(circleStr,cPathStr,I)

	'SetCoilWire(h,ra,ri,phi,N,O,wirePath)
	'ThickenWire (wirePath,c,cR)

	' older version w/o inputs
	'' SetCoilWireCurrent
	' SetCurrentPath
	' SetCoilWire
	' ThickenWire

End Sub

Sub SetUnit ()
	'set the units
	With Units
	    .Geometry "mm"
	    .Frequency "kHz"
	    .Voltage "V"
	    .Resistance "Ohm"
	    .Inductance "NanoH"
	    .TemperatureUnit  "Kelvin"
	    .Time "ns"
	    .Current "A"
	    .Conductance "Siemens"
	    .Capacitance "PikoF"
	End With

End Sub

' Sub SetSolver(F As String, Fnum As Integer)
Sub SetSolver(ByRef  F() As String,  Fstart As Double, Fend As Double, Fstep As Double, Fnum As Integer, cc As Integer)
	ChangeSolverType "LF Frequency Domain"
	With LFSolver
	    .ResetFrequencySettings
		For cc = 0 To Fnum-1
			F(cc) = cstr(Fstart + cc*Fstep)
			' MsgBox F(cc)
			.AddFrequency cstr(Fstart + cc*Fstep)
		Next cc

	End With

	With LFSolver
	     .Reset
	     .Method "Tetrahedral Mesh"
	     .Accuracy "1e-6"
	     .MaxLinIter "0"
	     .CalcImpedanceMatrix "False"
	     .StoreResultsInCache "False"
	     .Preconditioner "ILU"
	     .MeshAdaption "False"
	     .EquationType "Magnetoquasistatic"
	     .ValueScaling "rms"
	     .LSESolverType "Auto"
	     .SetTreeCotreeGauging "True"
	     .EnableDivergenceCheck "True"
	     .TetSolverOrder "2"
	     .TetAdaption "False"
	     .UseMaxNumberOfThreads "True"
	     .MaxNumberOfThreads "72"
	     .MaximumNumberOfCPUDevices "2"
	     .UseDistributedComputing "False"
	End With
End Sub

Sub SetBoundBox (xminb As Double,yminb As Double,zminb As Double,xmaxb As Double,ymaxb As Double,zmaxb As Double)
	Plot.DrawBox True
	With Background
	     .ResetBackground
	     .Type "Normal"
	     .Epsilon "1.0"
	     .Mue "1.0"
	     .Rho "1.204"
	     .ThermalType "Normal"
	     .ThermalConductivity "0.026"
	     .HeatCapacity "1.005"
	     .XminSpace xminb * Units.GetGeometrySIToUnit()
	     .XmaxSpace xmaxb * Units.GetGeometrySIToUnit()
	     .YminSpace yminb * Units.GetGeometrySIToUnit()
	     .YmaxSpace ymaxb * Units.GetGeometrySIToUnit()
	     .ZminSpace zminb * Units.GetGeometrySIToUnit()
	     .ZmaxSpace zmaxb * Units.GetGeometrySIToUnit()
	     .ApplyInAllDirections "True"
	End With
	With Boundary
	     .Xmin "electric"
	     .Xmax "electric"
	     .Ymin "electric"
	     .Ymax "electric"
	     .Zmin "electric"
	     .Zmax "electric"
	     .Xsymmetry "none"
	     .Ysymmetry "none"
	     .Zsymmetry "none"
	End With
	With Mesh
	     .MeshType "Tetrahedral"
	     .RatioLimit "50"
	     .MinimumStepNumber "20"
	     .MinimumStepNumberTet "5"
	     .SetAutomeshRefineDielectricsType "None"
	     .DensityTransitions "0.8"
	End With
	With MeshSettings
	     .SetMeshType "Hex"
	     .Set "StepsPerBoxNear", "20"
	     .Set "StepsPerBoxFar", "10"
	     .SetMeshType "Tet"
	     .Set "StepsPerBoxFar", "5"
	     .SetMeshType "Plane"
	     .Set "StepsPerBoxFar", "5"
	End With
	With MeshAdaption3D
	     .SetLFAdaptionStrategy "Energy"
	End With
End Sub

Sub SetCoilWire (h As String,ra As String,ri As String,phi As String,N As String,O As String,wirePath As String)
' Sub SetCoilWire (h As Double,ra As Double,ri As Double,phi As Double,N As Double,O As Double)

	Dim scst_torrus_ri As Double, scst_torrus_ra As Double,cst_xxx As Double, scst_torrus_phi As Double, scst_torrus_h As Double
	Dim scst_torrus_N As Integer, cst_result As Integer,  sCurveName As String, scst_clock As Integer, cst_clock As Integer
	Dim cst_torrus_N As Double, cst_torrus_phi As Double, cst_torrus_h As Double, cst_torrus_ra As Double, cst_torrus_ri As Double
	Dim first As String,  seconditem As String,cName As String, clock_yes_no As Integer, cClock As Integer
	cst_result = -1%
	' cName = "Linear_Spiral"
	cName = wirePath
	If (cst_result =0) Then Exit All   ' if cancel/help is clicked, exit all
	MakeSureParameterExists(cName+"_cst_torrus_h",h)
	MakeSureParameterExists(cName+"_cst_torrus_ra",ra)
	MakeSureParameterExists(cName+"_cst_torrus_ri",ri)
	MakeSureParameterExists(cName+"_cst_torrus_phi",phi)
	MakeSureParameterExists(cName+"_cst_torrus_N",N)
	MakeSureParameterExists(cName+"_cst_torrus_orientation",O)
	If (RestoreDoubleParameter(cName+"_cst_torrus_phi") <= 0) Then
		ReportWarningToWindow "Angle phi is not positive. Coil not constructed."
		Exit All
	End If
	 With Brick
			.Reset
			.Name "dummy_solid_trapez"
			.Component "Dummy_spiral_coil"
			.Material "PEC"
			.Xrange "0", cName+"_cst_torrus_h"
			.Yrange "0", cName+"_cst_torrus_ra"
			.Zrange "0", cName+"_cst_torrus_phi"
			.Create
		End With
		Component.Delete "Dummy_spiral_coil"
		cst_torrus_N= restoredoubleparameter (cName+"_cst_torrus_N")
		cst_torrus_h= restoredoubleparameter (cName+"_cst_torrus_h")
		cst_torrus_ra= restoredoubleparameter (cName+"_cst_torrus_ra")
		cst_torrus_ri= restoredoubleparameter (cName+"_cst_torrus_ri")
		cst_torrus_phi= restoredoubleparameter (cName+"_cst_torrus_phi")
		cst_clock = restoredoubleparameter (cName+"_cst_torrus_orientation")
	 ' Begin construction
	 clock_yes_no = CInt(cst_clock)
		On Error GoTo Curve_Exists
	 	Curve.NewCurve "3D-Linear-Spiral"
	 	Curve_Exists:
		On Error GoTo 0
	 	sCurveName = cName '"3dpolygon_1"
	 	With    Polygon3D
	  		.Reset
	  		.Name sCurveName
	  		.Curve "3D-Linear-Spiral"
	  		' .Curve "3D-Linear-Spiral"
	        ' the upper limit takes the numerical inaccuracies into account. The logic of the following loop is as follows:
	        ' We go complete the number of turns the user has specified. If 360 modulo cst_torrus_phi is not equal to zero
	        ' (which shouldn't be the case usually) we insert the last segment only if the overlap is less than half the length
	        ' of the segment.
	  		For cst_xxx = 0  To    cst_torrus_N *2*Pi+(cst_torrus_phi*pi/180)/2 STEP  cst_torrus_phi*pi/180
	  			If clock_yes_no = 1 Then
		 			.Point cst_torrus_ra*Sin(cst_xxx), cst_torrus_ri*Cos(cst_xxx) , cst_torrus_h*cst_xxx/(2*pi* cst_torrus_N)
		 		Else
		 			.Point -cst_torrus_ra*Sin(cst_xxx), cst_torrus_ri*Cos(cst_xxx) , cst_torrus_h*cst_xxx/(2*pi* cst_torrus_N)
		 		End If
	  		Next cst_xxx
	  		.Create
	 	End With
End Sub

Sub SetCurrent (circleStr As String, wirePath As String, I As String)
	With CurrentPath
	     .Reset
	     ' .Name "path1"
		 .Name circleStr
	     .Type "CurvePath"
	     '.Current "1"
		 .Current I
	     .Phase "0.0"
	     .PathCurve "3D-Linear-Spiral:" & wirePath
	     MsgBox "3D-Linear-Spiral:" & wirePath
	     .Add
	End With
End Sub



Sub SetCoilWireCurrent (h As String,ra As String,ri As String,phi As String,N As String,O As String, wirePath As String)

	Dim scst_torrus_ri As Double, scst_torrus_ra As Double,cst_xxx As Double, scst_torrus_phi As Double, scst_torrus_h As Double
	Dim scst_torrus_N As Integer, cst_result As Integer,  sCurveName As String, scst_clock As Integer, cst_clock As Integer
	Dim cst_torrus_N As Double, cst_torrus_phi As Double, cst_torrus_h As Double, cst_torrus_ra As Double, cst_torrus_ri As Double
	Dim first As String,  seconditem As String,cName As String, clock_yes_no As Integer, cClock As Integer
	cst_result = -1%
	' cName = "Linear_Spiral"
	cName = wirePath
	If (cst_result =0) Then Exit All   ' if cancel/help is clicked, exit all
	MakeSureParameterExists(cName+"_cst_torrus_h",h)
	MakeSureParameterExists(cName+"_cst_torrus_ra",ra)
	MakeSureParameterExists(cName+"_cst_torrus_ri",ri)
	MakeSureParameterExists(cName+"_cst_torrus_phi",phi)
	MakeSureParameterExists(cName+"_cst_torrus_N",N)
	MakeSureParameterExists(cName+"_cst_torrus_orientation",O)
	If (RestoreDoubleParameter(cName+"_cst_torrus_phi") <= 0) Then
		ReportWarningToWindow "Angle phi is not positive. Coil not constructed."
		Exit All
	End If
	 With Brick
			.Reset
			.Name "dummy_solid_trapez"
			.Component "Dummy_spiral_coil"
			.Material "PEC"
			.Xrange "0", cName+"_cst_torrus_h"
			.Yrange "0", cName+"_cst_torrus_ra"
			.Zrange "0", cName+"_cst_torrus_phi"
			.Create
		End With
		Component.Delete "Dummy_spiral_coil"
		cst_torrus_N= restoredoubleparameter (cName+"_cst_torrus_N")
		cst_torrus_h= restoredoubleparameter (cName+"_cst_torrus_h")
		cst_torrus_ra= restoredoubleparameter (cName+"_cst_torrus_ra")
		cst_torrus_ri= restoredoubleparameter (cName+"_cst_torrus_ri")
		cst_torrus_phi= restoredoubleparameter (cName+"_cst_torrus_phi")
		cst_clock = restoredoubleparameter (cName+"_cst_torrus_orientation")

	 ' Begin construction
		clock_yes_no = CInt(cst_clock)
		On Error GoTo Curve_Exists
	 	Curve.NewCurve "3D-Linear-Spiral"
	 	Curve_Exists:
		On Error GoTo 0
	 	sCurveName = cName '"3dpolygon_1"

		' MsgBox "" & cst_torrus_ra
		' MsgBox "" & cst_torrus_ri
		' MsgBox "" & cst_torrus_h
		' MsgBox "" & cst_torrus_N
		' MsgBox "" & cst_torrus_phi
		' MsgBox "" & cst_torrus_N *2*Pi+(cst_torrus_phi*pi/180)/2
		' MsgBox "" & cst_torrus_phi*pi/180

	 	With    Polygon3D
	  		.Reset
	  		.Name sCurveName
	  		.Curve "3D-Linear-Spiral"
	        ' the upper limit takes the numerical inaccuracies into account. The logic of the following loop is as follows:
	        ' We go complete the number of turns the user has specified. If 360 modulo cst_torrus_phi is not equal to zero
	        ' (which shouldn't be the case usually) we insert the last segment only if the overlap is less than half the length
	        ' of the segment.
			Dim helixFINAL As Double, helixSTEP As Double
			helixFINAL = cst_torrus_N *2*Pi+(cst_torrus_phi*pi/180)/2
			helixSTEP = cst_torrus_phi*pi/180
	  		For cst_xxx = 0  To   helixFINAL STEP  helixSTEP
	  			If clock_yes_no = 1 Then
		 			.Point cst_torrus_ra*Sin(cst_xxx), cst_torrus_ri*Cos(cst_xxx) , cst_torrus_h*cst_xxx/(2*pi* cst_torrus_N)
		 		Else
		 			.Point -cst_torrus_ra*Sin(cst_xxx), cst_torrus_ri*Cos(cst_xxx) , cst_torrus_h*cst_xxx/(2*pi* cst_torrus_N)
		 		End If
	  		Next cst_xxx
	  		.Create
	End With
End Sub


Sub ThickenWire (wirePath As String, c As String, cR As String)
	' Set
	' Pick.PickCurveEndpointFromId "3D-Linear-Spiral:Linear_Spiral", "1"
	Pick.PickCurveEndpointFromId wirePath, "1"
	WCS.AlignWCSWithSelectedPoint
	With WCS
	     .SetNormal "1", "0", "0"
	     .SetUVector "0", "1", "0"
	     .ActivateWCS "local"
	End With

	' this will have to dynamincally changed
	With Circle
	     .Reset
		 .Name c
	     '.Name "circle1"
	     .Curve "3D-Linear-Spiral"
	     ' .Radius "0.075"
		 .Radius cR
	     .Xcenter "0.005"
     	.Ycenter "-0"
     	.Segments "0"
    	.Create
	End With

	Component.New "component1"
	With SweepCurve
	     .Reset
	     .Name "solid1"
	     .Component "component1"
	     .Material "PEC"
	     .Twistangle "0.0"
	     .Taperangle "0.0"
	     .ProjectProfileToPathAdvanced "True"
	     .DeleteProfile "True"
	     .DeletePath "True"
 		 .Path "3D-Linear-Spiral:" & wirePath
	     ' .Path "3D-Linear-Spiral:Linear_Spiral"
		 .Curve "3D-Linear-Spiral:" & c
	     '.Curve "3D-Linear-Spiral:circle1"
	     .Create
	End With
End Sub

