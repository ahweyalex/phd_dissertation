'#Language "WWB-COM"

Option Explicit

Sub Main
	'1) initialize set up
	SetUnit
	' Set up Frequency Sweep
	Dim Fstart As Double, Fend As Double, Fstep As Double
	Dim Fnum As Integer, cc As Integer, F() As String
		Fend = 125
		Fstart = 125
		Fstep = 1
		Fnum = Int(((Fend - Fstart)/Fstep) + 1)
	ReDim F(Fnum)
	SetSolver(F,Fstart,Fend,Fstep,Fnum,cc)

	'2) Set current direction(path)
	' Create helix wire
	Dim h As String,ra As String,ri As String,phi As String,N As String,O As String, cPathStr As String, wPathStr As String, wT As String
		' wT = "0.00051"
		wT = "0.25"
		ra = "1"
		ri = "1"
		phi = "10"
		N = "4"
		O = "1"
		h  = 1.1*CStr(2*CDbl(wT)*CDbl(N))
		' MsgBox h
		cPathStr = "currentWire"
		wPathStr = "Wire"
	SetCoilWireCurrent(h,ra,ri,phi,N,O,cPathStr)
	Dim circleStr As String, I As String
	'Set current along helix wire
		circleStr = "c"
		I = "1"
	SetCurrent(circleStr,cPathStr,I)
	' 3) Set vacuum objects around current(path) to increase meshing
	SetCoilWireCurrent(h,ra,ri,phi,N,O,cPathStr)
	ThickenWire(cPathStr,wT)
	'4) Set up Boundary Box
	Dim xminb, yminb, zminb, xmaxb, ymaxb, zmaxb
		xminb = CDbl(h)
		yminb = CDbl(h)
		zminb = CDbl(h)
		xmaxb = CDbl(h)
		ymaxb = CDbl(h)
		zmaxb = CDbl(h)
	SetBoundBox(xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)


End Sub

Sub SetUnit ()
	'set the units
	With Units
	    .Geometry "m"
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
	UseDistributedComputingForParameters "False"
	MaxNumberOfDistributedComputingParameters "2"
	UseDistributedComputingMemorySetting "False"
	MinDistributedComputingMemoryLimit "0"
	UseDistributedComputingSharedDirectory "False"

End Sub

Sub SetBoundBox (xminb As Double,yminb As Double,zminb As Double,xmaxb As Double,ymaxb As Double,zmaxb As Double)
	With Background
	     .ResetBackground
	     .XminSpace xminb
	     .XmaxSpace xmaxb
	     .YminSpace yminb
	     .YmaxSpace ymaxb
	     .ZminSpace zminb
	     .ZmaxSpace zmaxb
	     .ApplyInAllDirections "True"
	End With

	With Material
	     .Reset
	     .Rho "1.204"
	     .ThermalType "Normal"
	     .ThermalConductivity "0.026"
	     .HeatCapacity "1.005"
	     .DynamicViscosity "1.84e-5"
	     .Emissivity "0.0"
	     .MetabolicRate "0.0"
	     .VoxelConvection "0.0"
	     .BloodFlow "0"
	     .MechanicsType "Unused"
	     .FrqType "all"
	     .Type "Normal"
	     .MaterialUnit "Frequency", "Hz"
	     .MaterialUnit "Geometry", "m"
	     .MaterialUnit "Time", "s"
	     .MaterialUnit "Temperature", "Kelvin"
	     .Epsilon "1.00059"
	     .Mu "1.0"
	     .Sigma "0.0"
	     .TanD "0.0"
	     .TanDFreq "0.0"
	     .TanDGiven "False"
	     .TanDModel "ConstSigma"
	     .EnableUserConstTanDModelOrderEps "False"
	     .ConstTanDModelOrderEps "1"
	     .SetElParametricConductivity "False"
	     .ReferenceCoordSystem "Global"
	     .CoordSystemType "Cartesian"
	     .SigmaM "0"
	     .TanDM "0.0"
	     .TanDMFreq "0.0"
	     .TanDMGiven "False"
	     .TanDMModel "ConstSigma"
	     .EnableUserConstTanDModelOrderMu "False"
	     .ConstTanDModelOrderMu "1"
	     .SetMagParametricConductivity "False"
	     .DispModelEps  "None"
	     .DispModelMu "None"
	     .DispersiveFittingSchemeEps "Nth Order"
	     .MaximalOrderNthModelFitEps "10"
	     .ErrorLimitNthModelFitEps "0.1"
	     .UseOnlyDataInSimFreqRangeNthModelEps "False"
	     .DispersiveFittingSchemeMu "Nth Order"
	     .MaximalOrderNthModelFitMu "10"
	     .ErrorLimitNthModelFitMu "0.1"
	     .UseOnlyDataInSimFreqRangeNthModelMu "False"
	     .UseGeneralDispersionEps "False"
	     .UseGeneralDispersionMu "False"
	     .NLAnisotropy "False"
	     .NLAStackingFactor "1"
	     .NLADirectionX "1"
	     .NLADirectionY "0"
	     .NLADirectionZ "0"
	     .Colour "0.6", "0.6", "0.6"
	     .Wireframe "False"
	     .Reflection "False"
	     .Allowoutline "True"
	     .Transparentoutline "False"
	     .Transparency "0"
	     .ChangeBackgroundMaterial
	End With

	With Boundary
     .Xmin "open"
     .Xmax "open"
     .Ymin "open"
     .Ymax "open"
     .Zmin "open"
     .Zmax "open"
     .Xsymmetry "none"
     .Ysymmetry "none"
     .Zsymmetry "none"
     .ApplyInAllDirections "True"
	End With

	With Mesh
	     .MeshType "Tetrahedral"
	     .SetCreator "Low Frequency"
	End With
	With MeshSettings
	     .SetMeshType "Tet"
	     .Set "Version", 1%
	     'MAX CELL - WAVELENGTH REFINEMENT
	     .Set "StepsPerWaveNear", "4"
	     .Set "StepsPerWaveFar", "4"
	     .Set "PhaseErrorNear", "0.02"
	     .Set "PhaseErrorFar", "0.02"
	     .Set "CellsPerWavelengthPolicy", "automatic"
	     'MAX CELL - GEOMETRY REFINEMENT
	     .Set "StepsPerBoxNear", "20"
	     .Set "StepsPerBoxFar", "5"
	     .Set "ModelBoxDescrNear", "maxedge"
	     .Set "ModelBoxDescrFar", "maxedge"
	     'MIN CELL
	     .Set "UseRatioLimit", "0"
	     .Set "RatioLimit", "100"
	     .Set "MinStep", "0"
	     'MESHING METHOD
	     .SetMeshType "Unstr"
	     .Set "Method", "0"
	End With
	With MeshSettings
	     .SetMeshType "Tet"
	     .Set "CurvatureOrder", "1"
	     .Set "CurvatureOrderPolicy", "automatic"
	     .Set "CurvRefinementControl", "NormalTolerance"
	     .Set "NormalTolerance", "22.5"
	     .Set "SrfMeshGradation", "1.5"
	     .Set "SrfMeshOptimization", "1"
	End With
	With MeshSettings
	     .SetMeshType "Unstr"
	     .Set "UseMaterials",  "0"
	     .Set "MoveMesh", "0"
	End With
	With MeshSettings
	     .SetMeshType "Tet"
	     .Set "UseAnisoCurveRefinement", "1"
	     .Set "UseSameSrfAndVolMeshGradation", "1"
	     .Set "VolMeshGradation", "1.5"
	     .Set "VolMeshOptimization", "1"
	End With
	With MeshSettings
	     .SetMeshType "Unstr"
	     .Set "SmallFeatureSize", "0"
	     .Set "CoincidenceTolerance", "1e-06"
	     .Set "SelfIntersectionCheck", "1"
	     .Set "OptimizeForPlanarStructures", "0"
	End With
	With Mesh
	     .SetParallelMesherMode "Tet", "maximum"
	     .SetMaxParallelMesherThreads "Tet", "1"
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
	     ' MsgBox "3D-Linear-Spiral:" & wirePath
	     .Add
	End With
End Sub

Sub SetCoilWireCurrent (h As String,ra As String,ri As String,phi As String,N As String,O As String, wirePath As String)
	' declare variables
	Dim scst_torrus_ri As Double, scst_torrus_ra As Double,cst_xxx As Double, scst_torrus_phi As Double, scst_torrus_h As Double
	Dim scst_torrus_N As Integer, cst_result As Integer,  sCurveName As String, scst_clock As Integer, cst_clock As Integer
	Dim cst_torrus_N As Double, cst_torrus_phi As Double, cst_torrus_h As Double, cst_torrus_ra As Double, cst_torrus_ri As Double
	Dim first As String,  seconditem As String,cName As String, clock_yes_no As Integer, cClock As Integer

	cst_result = -1%
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

Sub ThickenWire(wireName As String,wT As String)
	With Wire
	     .Reset
	     .Name "wire1"
	     .Folder "3D-Linear-Spiral"
	     .Radius wT
	     .Type "CurveWire"
	     .Curve "3D-Linear-Spiral:" & wireName
		 .Material "Copper (annealed)"
	     ' .Material "Vacuum"
	     .SolidWireModel "True"
	     .Termination "Natural"
	     .Mitering "NewMiter"
	     .AdvancedChainSelection "True"
	     .Add
	End With

	With Wire
     .Reset
     .SolidName "component1:wire1"
     .Name "wire1"
     .Folder "3D-Linear-Spiral"
     .KeepWire "False"
     .ConvertToSolidShape
	End With

End Sub
