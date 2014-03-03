# -*- coding: utf-8 -*-
# Density Map Module of VEDA
#
# Copyright:
# 2009-01-01 Gael Goret  gael.goret@ibs.fr
#
# Last modified:
# 2010-09-06 Gael Goret  gael.goret@ibs.fr
import vtk
try :
	import Tkinter as tk
except :
	import tkinter as tk
from math import sqrt
import gfx,Map,mod,itf,sym
from os import system,chdir,path
import tkColorChooser
import tkMessageBox as MB
#===============================================================================
class Map():
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def __init__(self):
		self.id=None
		self.fn=None
		self.acteur=None
		self.reader=None
		self.mapper=None
		self.iso=None
		self.sigma=1.0
		self.avg=0.
		self.isov=1.0
		self.box=None
		self.rendtype='wireframe'
		self.color=(0,0.5,0.75)
		self.scale=1.0
		self.oldscale=1.0
		self.opct=1.0
		self.isdeci='0'
		self.issmooth='1'
		self.ratio=1.0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def rescale_map(self,gfx,scale,caller):
		cropentry=None
		nfv=None
		load_map(gfx,0,gfx.itf.root,gfx.itf.status,scale,self.rendtype,self.isov,self.opct,cropentry,nfv,self.isdeci,self.issmooth,self.color,caller)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def map_res_range(self):
		minres = max(self.reader.GetOutput().GetSpacing())*2+0.05
		(xmin,xmax,ymin,ymax,zmin,zmax)=self.box.GetBounds()
		maxres = min(xmax-xmin,ymax-ymin,zmax-zmin)
		return (maxres,minres)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def assign_map(gfx,nfv,cpltnfv,scaleentry,clentry,opctentry,butcolor):
	fn=itf.browsemap()
	if fn =='' or fn == None or fn ==() :
		if gfx.map!=[]:
			nfv.set(extract_file_from_path(gfx.map[0].fn))
			cpltnfv.set(gfx.map[0].fn)
		else :
			nfv.set('')
			cpltnfv.set('')
	else :
		cpltnfv.set(fn)
		nfv.set(extract_file_from_path(fn))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def extract_file_from_path(fn):
	for i in range(len(fn)):
		if fn[-(i+1)]=='/':
			return fn[-i:]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def color_map(gfx,butcolor):
	'''
	if gfx.map==[]:
		MB.showwarning('Info','Select map file')
		return
	'''
	if gfx.map!=[]:
		new_color = tkColorChooser.askcolor (title="Map color",initialcolor=mod.vtk2tk_color(gfx.map[0].color))
	else:
		new_color = tkColorChooser.askcolor (title="Map color",initialcolor=mod.vtk2tk_color((0,0.5,0.75)))
	if new_color[1] != None:
		    	col = mod.tk2vtk_color (new_color[0])
    			butcolor.config(bg=mod.vtk2tkhex_color(col))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def set_map_id(gfx):
	if gfx.map==[]:
		return 0
	else:
		return gfx.map[-1].id + 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def init_cam_slab(gfx,bounds):
	""" set initial slab using dimension of box around EM map """
	(xmin,xmax,ymin,ymax,zmin,zmax)=bounds
	maxi=max(xmax-xmin,ymax-ymin,zmax-zmin)
	mini=min(xmax-xmin,ymax-ymin,zmax-zmin)
	fp=((xmax+xmin)/2.,(ymax+ymin)/2.,(zmax+zmin)/2.)
	gfx.camera.SetFocalPoint(fp[0],fp[1],fp[2])
	gfx.camera.SetPosition(fp[0],fp[1],fp[2]+(zmax-zmin)*4.)
	s=gfx.slab=zmax-zmin
	z=gfx.zpos=0
	c=gfx.camera.GetDistance()
	sm=gfx.slabmax=2*s
	zm=gfx.zposmax=2.2*s
	gfx.itf.init_slab_stuff(gfx,s*1.5,z,sm,zm)
	#gfx.itf.init_fog_stuff(gfx,s*1.5,z,sm,zm)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Chargement de la map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_map(gfx,mapfile,root,status,scale,rendtype,isov,opct,cropentry,nfv,vardeci,varsmooth,color,caller):
	if caller != 'fit':
		root.configure(cursor='watch')
		status.set('Map loading ... please wait')
	if mapfile =='' or mapfile == None or mapfile ==() :
		MB.showwarning('Info','Select map file')
		status.clear()
		root.configure(cursor='arrow')
		return
	try:
		gfx.renderer.RemoveActor(gfx.map[0].acteur)
		gfx.renderer.RemoveActor(gfx.map[0].box)
	except:
		pass
	if mapfile == 0:
		mapfile = gfx.map[0].fn
	if gfx.map == []:
		gfx.map = [Map()]
		gfx.map[0].id=0
	if 'map' in caller :
		clean_map(gfx) #supression des fichiers sort.s xudi et iudi
	if '.vtk' in mapfile:
		chdir(gfx.tmpdir)
		v2v_out='info_map'
		if caller != 'crop':
			gfx.map[0].sigma,gfx.map[0].avg = map_sigma_avg(v2v_out)
		if scale != gfx.map[0].scale:
			spc = None
			o = None
			f = open(mapfile,'r')
			for l in f:
				if l.startswith('SPACING'):
					spc = l.split()[1:4]
				if l.startswith('ORIGIN'):
					o = l.split()[1:4]
				if spc != None and o != None:
					break
			f.close()
			gfx.map[0].ratio = scale/gfx.map[0].scale
			if spc != None and o != None:
				system("sed -i -e /^SPACING/s/.*/'SPACING %f %f %f'/ %s"%(float(spc[0])*gfx.map[0].ratio,float(spc[1])*gfx.map[0].ratio,float(spc[2])*gfx.map[0].ratio,mapfile))
				system("sed -i -e /^ORIGIN/s/.*/'ORIGIN %f %f %f'/ %s"%(float(o[0])*gfx.map[0].ratio,float(o[1])*gfx.map[0].ratio,float(o[2])*gfx.map[0].ratio,mapfile))
		chdir(gfx.workdir)
	if '.ezd' in mapfile:
		chdir(gfx.tmpdir)
		mapfileout = extract_file_from_path(mapfile)[:-4]+'.vtk'
		e2v_out='info_map'
		system(gfx.vedabin+'/e2v.exe >> %s <<ENDOF\n%s  \n%f  \n%s  \nENDOF'%(e2v_out,mapfile,scale,mapfileout))
		mapfile = gfx.tmpdir + '/' + mapfileout
		gfx.map[0].sigma,gfx.map[0].avg = map_sigma_avg(e2v_out)
		chdir(gfx.workdir)
	gfx.map[0].fn = mapfile
	gfx.map[0].id = set_map_id(gfx)
	gfx.map[0].color = color
	gfx.map[0].oldscale = gfx.map[0].scale
	gfx.map[0].scale = scale
	if nfv !=None:
		nfv.set(extract_file_from_path(gfx.map[0].fn))
	reader = vtk.vtkStructuredPointsReader()
	reader.SetFileName(mapfile)
	reader.Update() #by calling Update() we read the file
	gfx.map[0].reader=reader
	iso = vtk.vtkMarchingContourFilter()
	iso.UseScalarTreeOn()
	iso.ComputeNormalsOn()
	iso.SetInputConnection(reader.GetOutputPort())
	iso.SetValue(0,isov*gfx.map[0].sigma+gfx.map[0].avg)
	gfx.map[0].iso=iso
	gfx.map[0].isov=isov
	if varsmooth == '1':
		#generate vectors
		clean = vtk.vtkCleanPolyData()
	  	clean.SetInputConnection(iso.GetOutputPort())
	 	clean.ConvertStripsToPolysOn()
		smooth = vtk.vtkWindowedSincPolyDataFilter()
	 	smooth.SetInputConnection(clean.GetOutputPort())
	 	smooth.BoundarySmoothingOn()
	 	smooth.GenerateErrorVectorsOn()
	  	smooth.GenerateErrorScalarsOn()
	  	smooth.NormalizeCoordinatesOn()
	  	smooth.NonManifoldSmoothingOn()
	  	smooth.FeatureEdgeSmoothingOn()
	  	smooth.SetEdgeAngle(90)
		smooth.SetFeatureAngle(90)
		smooth.Update()
	if vardeci=='1':
		deci = vtk.vtkDecimatePro()
		if varsmooth == '0':
			deci.SetInput(iso.GetOutput())
		else :
			deci.SetInput(smooth.GetOutput())
		deci.PreserveTopologyOn()
		deci.BoundaryVertexDeletionOn()
		deci.SplittingOn()
		deci.PreSplitMeshOn()
		deci.SetTargetReduction(0.97)
		gfx.map[0].isdeci='1'
		mapper = vtk.vtkOpenGLPolyDataMapper()
		mapper.SetInputConnection(deci.GetOutputPort())
	else :
		mapper = vtk.vtkOpenGLPolyDataMapper()
		if varsmooth == '1':
			mapper.SetInputConnection(smooth.GetOutputPort()) ### <- connection here
		else :
			mapper.SetInputConnection(iso.GetOutputPort())
		#mapper.SetInput(newpd) ### <- newpd connect there
		gfx.map[0].isdeci='0'
	mapper.ScalarVisibilityOff()
	mapper.Update()
	gfx.map[0].mapper=mapper
	actor = vtk.vtkOpenGLActor()
	actor.SetMapper(mapper)
	gfx.map[0].acteur=actor
	#actor.SetScale(scale,scale,scale) gerer differament
	actor.GetProperty().SetColor(gfx.map[0].color)
	actor.PickableOff()
	#definition de la box
	outline = vtk.vtkOutlineFilter()
        outline.SetInput(reader.GetOutput())
        outlineMapper = vtk.vtkPolyDataMapper()
        outlineMapper.SetInput(outline.GetOutput())
	box=vtk.vtkActor()
        box.SetMapper( outlineMapper )
        box.GetProperty().SetColor((invcolor(gfx.map[0].color)))
        box.PickableOff()
	#box.SetScale(scale,scale,scale)
	gfx.map[0].box = box
	#get boxwidget bounds and set axes lenth
	(xmin,xmax,ymin,ymax,zmin,zmax)=box.GetBounds()
        x=abs(xmin-xmax)/2.0
	y=abs(ymin-ymax)/2.0
	z=abs(zmin-zmax)/2.0
	gfx.axes.SetTotalLength( x, y , z ) #defini la longeurs des axe
	init_cam_slab(gfx,(xmin,xmax,ymin,ymax,zmin,zmax)) #defini le slab correct
	gfx.map[0].rendtype=rendtype
	if rendtype=='Wireframe':
		actor.GetProperty().SetRepresentationToWireframe()
	elif rendtype=='Surface':
		actor.GetProperty().SetRepresentationToSurface()
	elif rendtype=='Points':
		actor.GetProperty().SetRepresentationToPoints()
		actor.GetProperty().SetPointSize(5)
	else :
		actor.GetProperty().SetRepresentationToWireframe()
	gfx.map[0].opct=opct
	actor.GetProperty().SetOpacity(opct)
	actor.GetProperty().SetInterpolationToGouraud()
	actor.GetProperty().SetSpecular(.4)
	actor.GetProperty().SetSpecularPower(10)
	if cropentry!=None:
		if gfx.crop==None:
			gfx.crop=Crop(gfx,iso,cropentry,None) #here entryval = None
	rendermap(gfx)
	#ajustement pour la symetry helicoidale
	if gfx.map[0].scale != gfx.map[0].oldscale:#changement de scale
		gfx.itf.propagate_scale(gfx,scale,caller)
		if gfx.ps != None:
			if gfx.ps.solidtype == 'Helicoidal':
				gfx.ps.display_tube(gfx,caller)
			elif gfx.ps.solidtype == 'Icosahedral' or gfx.ps.solidtype =='Octahedral' or gfx.ps.solidtype == 'Tetrahedral':
				gfx.ps.display_platonic(gfx,gfx.ps.ori)
			elif gfx.ps.solidtype =='Cn' or gfx.ps.solidtype == 'Dn':
				gfx.ps.display_Xn(gfx)
	if caller == 'crop':#crop uniquement helicoidal
		if gfx.ps != None:
			if gfx.ps.solidtype == 'Helicoidal':
				gfx.ps.display_tube(gfx,caller)
	if caller != 'fit':
		status.clear()
		root.configure(cursor='arrow')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def clean_map(gfx):
	chdir(gfx.workdir)
	if path.exists('d/data.d'):
		system('rm -f d/data.d')
	if path.exists('d/emap.d'):
		system('rm -f d/emap.d')
	if path.exists('o/sort.s'):
		system('rm -f o/sort.s')
	if path.exists('f/xudi'):
		system('rm -f f/xudi')
	if path.exists('f/iudi'):
		system('rm -f f/iudi')
	if gfx.map != []:
		if gfx.map[0].fn!=None:
			if path.exists(gfx.map[0].fn):
				system('rm -f %s'%gfx.map[0].fn)
	gfx.crop = None
	gfx.fit = None
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def map_sigma_avg(info_map):
	try:
		f=open(info_map,"r")
	except IOError:
		return (1.,0.)
	for l in f:
		if 'sigma' in l:
			sig, avg = float(l.split()[-1]),float(l.split()[-2])
	f.close()
	return (sig,avg)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              display de la map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def display_map(gfx,rendtype,isov,opct,color):
	if gfx.map==[]:
		MB.showwarning('Info','Load a map')
		return
	gfx.itf.root.configure(cursor='watch')
	gfx.itf.status.set('Rendering in progress ... please wait')
	gfx.map[0].iso.SetValue(0,float(isov)*gfx.map[0].sigma+gfx.map[0].avg)
	gfx.map[0].mapper.Update()
	gfx.map[0].color = color
	gfx.map[0].acteur.GetProperty().SetColor(gfx.map[0].color)
        gfx.map[0].box.GetProperty().SetColor((invcolor(gfx.map[0].color)))
  	gfx.map[0].rendtype=rendtype
	if rendtype=='Wireframe':
		gfx.map[0].acteur.GetProperty().SetRepresentationToWireframe()
	elif rendtype=='Surface':
		gfx.map[0].acteur.GetProperty().SetRepresentationToSurface()
	elif rendtype=='Points':
		gfx.map[0].acteur.GetProperty().SetRepresentationToPoints()
		gfx.map[0].acteur.GetProperty().SetPointSize(5)
	else :
		actor.GetProperty().SetRepresentationToWireframe()
	gfx.map[0].opct=float(opct)
	gfx.map[0].acteur.GetProperty().SetOpacity(float(opct))
	gfx.map[0].acteur.GetProperty().SetInterpolationToGouraud()
	rendermap(gfx)
	gfx.itf.status.clear()
	gfx.itf.root.configure(cursor='arrow')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Rendu de la map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def rendermap(gfx):
	for Map in gfx.map:
		gfx.renderer.AddActor(Map.box)
		gfx.renderer.AddActor(Map.acteur)
	gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def hidemap(gfx):
	if gfx.map[0].acteur.GetVisibility():
		gfx.map[0].acteur.VisibilityOff()
	else:
		gfx.map[0].acteur.VisibilityOn()
	gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def hidemapbox(gfx):
	if gfx.map[0].box.GetVisibility():
		gfx.map[0].box.VisibilityOff()
	else:
		gfx.map[0].box.VisibilityOn()
	gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def invcolor(col):
	return (abs(col[0]-0.27),abs(col[1]-0.27),abs(col[2]-0.27))
'''
def invcolor(col):
	return (col[0],col[1],col[2])
'''
#===============================================================================
#                              Crop
#===============================================================================
class Crop(vtk.vtkBoxWidget):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def __init__(self,gfx,iso,cropentry,entryval):
		self.gfx = gfx
		self.entry = cropentry
		self.entryval = entryval
		self.crpfn=None
		self.SetProp3D(gfx.map[0].box)
		self.planes=vtk.vtkPlanes()
		self.SetInteractor(gfx.iren)
		self.SetInput(gfx.map[0].iso.GetOutput())
		self.SetPlaceFactor(1)
		self.PlaceWidget()
		self.InsideOutOn()
		self.SetRotationEnabled(0)
		self.GetPlanes(self.planes)
		self.AddObserver("EndInteractionEvent",self.SelectPolygons)
		self.AddObserver("InteractionEvent",self.Update_crop_bounds)
		self.inorout=1
		self.init_entrys()
		(xmin,xmax,ymin,ymax,zmin,zmax)=self.gfx.map[0].box.GetBounds()
		spcing=self.gfx.map[0].reader.GetOutput().GetSpacing()
		self.old=(nint(xmin/spcing[0]),nint(xmax/spcing[0]),nint(ymin/spcing[1]),nint(ymax/spcing[1]),nint(zmin/spcing[2]),nint(zmax/spcing[2]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def init_entrys(self):
		if self.entryval == None:
			(xmin,xmax,ymin,ymax,zmin,zmax)=self.gfx.map[0].box.GetBounds()
			spcing=self.gfx.map[0].reader.GetOutput().GetSpacing()
			(pxmin,pxmax,pymin,pymax,pzmin,pzmax)=(nint(xmin/spcing[0]),nint(xmax/spcing[0]),nint(ymin/spcing[1]),nint(ymax/spcing[1]),nint(zmin/spcing[2]),nint(zmax/spcing[2]))
			self.old=(pxmin,pxmax,pymin,pymax,pzmin,pzmax)
			self.entry[0].delete(0,tk.END)
			self.entry[0].insert(tk.END,str(pxmin))
			self.entry[1].delete(0,tk.END)
			self.entry[1].insert(tk.END,str(pxmax))
			self.entry[2].delete(0,tk.END)
			self.entry[2].insert(tk.END,str(pymin))
			self.entry[3].delete(0,tk.END)
			self.entry[3].insert(tk.END,str(pymax))
			self.entry[4].delete(0,tk.END)
			self.entry[4].insert(tk.END,str(pzmin))
			self.entry[5].delete(0,tk.END)
			self.entry[5].insert(tk.END,str(pzmax))
			self.entryval = [self.entry[i].get() for i in range(6)]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def reset_clipping(self):
		if self.gfx.map[0].mapper!=None:
			self.gfx.map[0].mapper.RemoveAllClippingPlanes()
		if self.gfx.mol != []:
			for mol in self.gfx.mol:
				mol.mapper.RemoveAllClippingPlanes()
		if self.gfx.ps != None:
			if self.gfx.ps.mapper!=None:
				self.gfx.ps.mapper.RemoveAllClippingPlanes()
			if self.gfx.ps.spheremapper!=None:
				self.gfx.ps.spheremapper.RemoveAllClippingPlanes()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def reset_crop_bounds(self):
		try :
			self.PlaceWidget(self.gfx.map[0].box.GetBounds())
		except AttributeError:
			return
		self.SelectPolygons('widget', 'event')
		self.gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def SelectPolygons(self,widget,event):
   		(bxmin,bxmax,bymin,bymax,bzmin,bzmax)=self.gfx.map[0].box.GetBounds()
		pd = vtk.vtkPolyData()
		self.GetPolyData(pd)
		(xmin,xmax,ymin,ymax,zmin,zmax)=pd.GetBounds()
		#(xo,yo,zo)= pd.GetCenter()
		tmpxmin=max(bxmin,xmin)
		tmpxmax=min(bxmax,xmax)
		if tmpxmin < tmpxmax :
			xmin=tmpxmin
			xmax=tmpxmax
		elif  tmpxmax <= bxmin: #a gauche
			xmax = bxmin+xmax-xmin
			xmin = bxmin
		elif  tmpxmin >= bxmax: #a droite
			xmin = bxmax-xmax+xmin
			xmax = bxmax
		else : print 'strange case xmin : %s | xmax : %s'%(xmin,xmax)
		tmpymin=max(bymin,ymin)
		tmpymax=min(bymax,ymax)
		if tmpymin < tmpymax :
			ymin=tmpymin
			ymax=tmpymax
		elif  tmpymax <= bymin: #a gauche
			ymax = bymin+ymax-ymin
			ymin = bymin
		elif  tmpymin >= bymax: #a droite
			ymin = bymax-ymax+ymin
			ymax = bymax
		else : print "strange case tudutudutudutu"
		tmpzmin=max(bzmin,zmin)
		tmpzmax=min(bzmax,zmax)
		if tmpzmin < tmpzmax :
			zmin=tmpzmin
			zmax=tmpzmax
		elif  tmpzmax <= bzmin: #a gauche
			zmax = bzmin+zmax-zmin
			zmin = bzmin
		elif  tmpzmin >= bzmax: #a droite
			zmin = bzmax-zmax+zmin
			zmax = bzmax
		else : print "strange case tudutudutudutu"
		spcing=self.gfx.map[0].reader.GetOutput().GetSpacing()
		(pxmin,pxmax,pymin,pymax,pzmin,pzmax)=(nint(xmin/spcing[0]),nint(xmax/spcing[0]),nint(ymin/spcing[1]),nint(ymax/spcing[1]),nint(zmin/spcing[2]),nint(zmax/spcing[2])) #transforme en unite pixel
		##############################################
		(bxl,bxu,byl,byu,bzl,bzu) = self.gfx.map[0].box.GetBounds()
		(bxl,bxu,byl,byu,bzl,bzu) = (nint(bxl/spcing[0]),nint(bxu/spcing[0]),nint(byl/spcing[1]),nint(byu/spcing[1]),nint(bzl/spcing[2]),nint(bzu/spcing[2]))
		extx,exty,extz = pxmax-pxmin+1,pymax-pymin+1,pzmax-pzmin+1
		(xlow,xup) = (pxmin,pxmax)
		mextx = self.Control_extent(extx,'-')
		deltam=abs(mextx-extx)
		pextx = self.Control_extent(extx,'+')
		deltap=abs(pextx-extx)
		if pxmin != self.old[0]:
			xlow = pxmin + (extx - pextx)
		if pxmax != self.old[1]:
			xup = pxmax - (extx - pextx)
		if (xlow < bxl) or (xup > bxu) or (deltap > deltam) :
			if pxmin != self.old[0]:
				xlow = pxmin + (extx - mextx)
			if pxmax != self.old[1]:
				xup = pxmax - (extx - mextx)
		(pxmin,pxmax) = (xlow,xup)
		(ylow,yup) = (pymin,pymax)
		mexty = self.Control_extent(exty,'-')
		deltam=abs(mexty-exty)
		pexty = self.Control_extent(exty,'+')
		deltap=abs(pexty-exty)
		if pymin != self.old[2]:
			ylow = pymin + (exty - pexty)
		if pymax != self.old[3]:
			yup = pymax - (exty - pexty)
		if (ylow < byl) or (yup > byu) or (deltap > deltam) :
			if pymin != self.old[2]:
				ylow = pymin + (exty - mexty)
			if pymax != self.old[3]:
				yup = pymax - (exty - mexty)
		(pymin,pymax) = (ylow,yup)
		(zlow,zup) = (pzmin,pzmax)
		mextz = self.Control_extent(extz,'-')
		deltam=abs(mextz-extz)
		pextz = self.Control_extent(extz,'+')
		deltap=abs(pextz-extz)
		if pzmin != self.old[4]:
			zlow = pzmin + (extz - pextz)
		if pzmax != self.old[5]:
			zup = pzmax - (extz - pextz)
		if (zlow < bzl) or (zup > bzu) or (deltap > deltam) :
			if pzmin != self.old[4]:
				zlow = pzmin + (extz - mextz)
			if pzmax != self.old[5]:
				zup = pzmax - (extz - mextz)
		(pzmin,pzmax) = (zlow,zup)
		self.old=(pxmin,pxmax,pymin,pymax,pzmin,pzmax)
		##############################################
		self.entry[0].delete(0,tk.END)
		self.entry[0].insert(tk.END,str(pxmin))
		self.entry[1].delete(0,tk.END)
		self.entry[1].insert(tk.END,str(pxmax))
		self.entry[2].delete(0,tk.END)
		self.entry[2].insert(tk.END,str(pymin))
		self.entry[3].delete(0,tk.END)
		self.entry[3].insert(tk.END,str(pymax))
		self.entry[4].delete(0,tk.END)
		self.entry[4].insert(tk.END,str(pzmin))
		self.entry[5].delete(0,tk.END)
		self.entry[5].insert(tk.END,str(pzmax))
		self.entryval = [self.entry[i].get() for i in range(6)]
		(txmin,txmax,tymin,tymax,tzmin,tzmax)=(spcing[0]*pxmin,spcing[0]*pxmax,spcing[1]*pymin,spcing[1]*pymax,spcing[2]*pzmin,spcing[2]*pzmax)
		self.PlaceWidget(txmin,txmax,tymin,tymax,tzmin,tzmax)
		self.GetPlanes(self.planes)
   		self.gfx.map[0].mapper.SetClippingPlanes(self.planes)
		if self.gfx.mol != []:
			for mol in self.gfx.mol:
				mol.mapper.SetClippingPlanes(self.planes)
		'''
		if self.gfx.ps!=None:
			if self.gfx.ps.mapper!=None:
				self.gfx.ps.mapper.SetClippingPlanes(self.planes)
			if self.gfx.ps.spheremapper!=None:
				self.gfx.ps.spheremapper.SetClippingPlanes(self.planes)
		'''
		self.gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def Control_extent(self,n,signe):
		bnl=[ 46,  58,  62,  74,  82,  86,  92,  94, 106, 116, 118, 122,  \
		124, 134, 138, 142, 146, 148, 158, 164, 166, 172, 174, 178,  \
		184, 186, 188, 194, 202, 206, 212, 214, 218, 222, 226, 230,  \
		232, 236, 244, 246, 248, 254, 258, 262, 268, 274, 276, 278,  \
		282, 284, 290, 292, 296, 298, 302, 310, 314, 316, 318, 322,  \
		326, 328, 332, 334, 344, 346, 348, 354, 356, 358, 362, 366,  \
		368, 370, 372, 376, 382, 386, 388, 394, 398, 402, 404, 406,  \
		410, 412, 414, 422, 424, 426, 428, 430, 434, 436, 438, 444,  \
		446, 452, 454, 458, 460, 464, 466, 470, 472, 474, 478, 482,  \
		488, 492, 496, 498, 502, 506, 508, 514, 516, 518, 522, 524,  \
		526, 530, 534, 536, 538, 542, 548, 552, 554, 556, 558, 562,  \
		564, 566, 568, 574, 580, 582, 584, 586, 590, 592, 596, 598,  \
		602, 604, 606, 610, 614, 618, 620, 622, 626, 628, 632, 634,  \
		636, 638, 642, 644, 652, 654, 656, 658, 662, 664, 666, 668,  \
		670, 674, 678, 682, 688, 690, 692, 694, 696, 698, 706, 708,  \
		710, 712, 716, 718, 724, 730, 732, 734, 736, 738, 740, 742,  \
		744, 746, 752, 754, 758, 762, 764, 766, 772, 774, 776, 778,  \
		782, 786, 788, 790, 794, 796, 804, 806, 808, 812, 814, 820,  \
		822, 824, 826, 828, 830, 834, 844, 846, 848, 852, 854, 856,  \
		860, 868, 870, 872, 874, 876, 888, 890, 892, 894, 902, 904,  \
		906, 908, 916, 920, 928, 930, 932, 938, 940, 942, 944, 946,  \
		948, 954, 956, 962, 964, 966, 970, 976, 978, 984, 986, 992,  \
		994, 996]
		maxs=1000
		mcd=2
		n=(n/mcd)*mcd
		ndiv=0
		while n>maxs:
			n=n/mcd
			ndiv=ndiv+1
		n=(n/mcd)*mcd
		while n in bnl:
			if signe=='-':
				n=n-mcd
			if signe=='+':
				n=n+mcd
		if ndiv != 0 :
			n = n * mcd ** ndiv
		return n
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def Update_crop_bounds(self,widget,event):
		pd = vtk.vtkPolyData()
		self.GetPolyData(pd)
		(xmin,xmax,ymin,ymax,zmin,zmax)=pd.GetBounds()
		self.bounds=pd.GetBounds()
		(bxmin,bxmax,bymin,bymax,bzmin,bzmax)=self.gfx.map[0].box.GetBounds()
		spcing=self.gfx.map[0].reader.GetOutput().GetSpacing()
		(xmin,xmax,ymin,ymax,zmin,zmax)=(xmin/spcing[0],xmax/spcing[0],ymin/spcing[1],ymax/spcing[1],zmin/spcing[2],zmax/spcing[2])
		self.entry[0].delete(0,tk.END)
		self.entry[0].insert(tk.END,'%.2f'%xmin)
		self.entry[1].delete(0,tk.END)
		self.entry[1].insert(tk.END,'%.2f'%xmax)
		self.entry[2].delete(0,tk.END)
		self.entry[2].insert(tk.END,'%.2f'%ymin)
		self.entry[3].delete(0,tk.END)
		self.entry[3].insert(tk.END,'%.2f'%ymax)
		self.entry[4].delete(0,tk.END)
		self.entry[4].insert(tk.END,'%.2f'%zmin)
		self.entry[5].delete(0,tk.END)
		self.entry[5].insert(tk.END,'%.2f'%zmax)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def show(self):	#permet le rechargement du crop a la reouverture du menu map
		if self.entry[0].get() != "" and self.entry[1].get() != "" and self.entry[2].get() != "" and self.entry[3].get() != "" and self.entry[4].get()!="" and self.entry[5].get()!="" :
			(xmin,xmax,ymin,ymax,zmin,zmax)=(float(self.entry[0].get()),float(self.entry[1].get()),float(self.entry[2].get()),float(self.entry[3].get()),float(self.entry[4].get()),float(self.entry[5].get()))
			spcing=self.gfx.map[0].reader.GetOutput().GetSpacing()
			(txmin,txmax,tymin,tymax,tzmin,tzmax)=(spcing[0]*xmin,spcing[0]*xmax,spcing[1]*ymin,spcing[1]*ymax,spcing[2]*zmin,spcing[2]*zmax)
			self.PlaceWidget(txmin,txmax,tymin,tymax,tzmin,tzmax)
			self.GetPlanes(self.planes)
	   		self.gfx.map[0].mapper.SetClippingPlanes(self.planes)
			if self.gfx.mol != []:
				for mol in self.gfx.mol:
					mol.mapper.SetClippingPlanes(self.planes)
		self.SetEnabled(1)
		#self.gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def hide(self):
		self.SetEnabled(0)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def crop_file(self):
		spcing=self.gfx.map[0].reader.GetOutput().GetSpacing()
		(bxl,bxu,byl,byu,bzl,bzu) = self.gfx.map[0].box.GetBounds()
		(bxl,bxu,byl,byu,bzl,bzu) = (nint(bxl/spcing[0]),nint(bxu/spcing[0]),nint(byl/spcing[1]),nint(byu/spcing[1]),nint(bzl/spcing[2]),nint(bzu/spcing[2]))
		ent=[float(self.entry[i].get()) for i in range(6)]
		(pxmin,pxmax,pymin,pymax,pzmin,pzmax)=(nint(ent[0]),nint(ent[1]),nint(ent[2]),nint(ent[3]),nint(ent[4]),nint(ent[5]))
		extx,exty,extz = pxmax-pxmin+1,pymax-pymin+1,pzmax-pzmin+1
		(xlow,xup) = (pxmin,pxmax)
		mextx = self.Control_extent(extx,'-')
		deltam=abs(mextx-extx)
		pextx = self.Control_extent(extx,'+')
		deltap=abs(pextx-extx)
		if pxmin != bxl:
			xlow = pxmin + (extx - pextx)
		if pxmax != bxu:
			xup = pxmax - (extx - pextx)
		if (xlow < bxl) or (xup > bxu) or (deltap > deltam) :
			if pxmin != bxl:
				xlow = pxmin + (extx - mextx)
			if pxmax != bxu:
				xup = pxmax - (extx - mextx)
		(pxmin,pxmax) = (xlow,xup)
		(ylow,yup) = (pymin,pymax)
		mexty = self.Control_extent(exty,'-')
		deltam=abs(mexty-exty)
		pexty = self.Control_extent(exty,'+')
		deltap=abs(pexty-exty)
		if pymin != byl:
			ylow = pymin + (exty - pexty)
		if pymax != byu:
			yup = pymax - (exty - pexty)
		if (ylow < byl) or (yup > byu) or (deltap > deltam) :
			if pymin != byl:
				ylow = pymin + (exty - mexty)
			if pymax != byu:
				yup = pymax - (exty - mexty)
		(pymin,pymax) = (ylow,yup)
		(zlow,zup) = (pzmin,pzmax)
		mextz = self.Control_extent(extz,'-')
		deltam=abs(mextz-extz)
		pextz = self.Control_extent(extz,'+')
		deltap=abs(pextz-extz)
		if pzmin != bzl:
			zlow = pzmin + (extz - pextz)
		if pzmax != bzu:
			zup = pzmax - (extz - pextz)
		if (zlow < bzl) or (zup > bzu) or (deltap > deltam) :
			if pzmin != bzl:
				zlow = pzmin + (extz - mextz)
			if pzmax != bzu:
				zup = pzmax - (extz - mextz)
		(pzmin,pzmax) = (zlow,zup)
		self.entry[0].delete(0,tk.END)
		self.entry[0].insert(tk.END,'%d'%pxmin)
		self.entry[1].delete(0,tk.END)
		self.entry[1].insert(tk.END,'%d'%pxmax)
		self.entry[2].delete(0,tk.END)
		self.entry[2].insert(tk.END,'%d'%pymin)
		self.entry[3].delete(0,tk.END)
		self.entry[3].insert(tk.END,'%d'%pymax)
		self.entry[4].delete(0,tk.END)
		self.entry[4].insert(tk.END,'%d'%pzmin)
		self.entry[5].delete(0,tk.END)
		self.entry[5].insert(tk.END,'%d'%pzmax)
		(txmin,txmax,tymin,tymax,tzmin,tzmax)=(spcing[0]*pxmin,spcing[0]*pxmax,spcing[1]*pymin,spcing[1]*pymax,spcing[2]*pzmin,spcing[2]*pzmax)
		self.PlaceWidget(txmin,txmax,tymin,tymax,tzmin,tzmax)
		strbounds='%s %s %s %s %s %s'%(pxmin,pxmax,pymin,pymax,pzmin,pzmax)
		chdir(self.gfx.tmpdir)
		v2v_out='info_map'
		system(self.gfx.vedabin + '/v2v.exe >> %s <<ENDOF\n%s \n%s \n%s \nENDOF'%(v2v_out,self.gfx.map[0].fn,strbounds,self.gfx.map[0].fn[:-4]+'_Crp.vtk'))
		system('mv %s %s'%(self.gfx.map[0].fn[:-4]+'_Crp.vtk',self.gfx.map[0].fn))
		self.crpfn=self.gfx.map[0].fn
		chdir(self.gfx.workdir)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def nint(f):
	return int(round(f))
