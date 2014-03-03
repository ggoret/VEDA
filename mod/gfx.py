# -*- coding: utf-8 -*-
# Graphical Module of VEDA
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
	try :
		import tkinter as tk
	except :
		print 'VEDA : Error : Tkinter is not installed.'
import os,sys,time
import Map,mod,itf,Nma,uro,sym
from vtk.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor
import tkMessageBox as MB
import tkFileDialog as FD
try:
	import IPython
except ImportError:
	print 'Warning : Ipython not installed. \n VEDA will run anyway, but you will not get access to debug mode'
HAS_GL = 1
try:
	from OpenGL.GL import  glEnable,glDisable, glFogf, glFogi, GL_FOG, GL_FOG_MODE, GL_EXP2, glHint, GL_FOG_HINT, GL_NICEST, GL_LINEAR, GL_FOG_START, GL_FOG_END,GL_FOG_DENSITY
except ImportError:
	print 'Warning : Python-Opengl not installed. \n VEDA will run anyway, but you will not get depth cueing effect'
	HAS_GL = 0
#===============================================================================
class gfx():
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def __init__(self):
		#~~~~~~Graphism~~~~~~~~~~~~~~~~~~~~~~~~~#
		try:
			self.vedabin = os.environ['VEDA'] + "/" +os.environ['VEDA_BIN']
		except KeyError:
			print 'Error : Environment variables VEDA or/and VEDA_BIN not defined'
			return
		self.version = "22.12.10"
		self.vedadir = sys.path[0][:-3]
		self.workdir = os.getcwd()
		self.tmpdir = self.create_tmp_dir()
		self.itf = None
		self.ifit = None #interactive fit
		self.fit = None #refinement fitting
		self.ps = None #Platonic Solid which support symmetry
		self.crop = None #crop object (widgetbox + vtkplanes)
		self.sym = None #constellation objet
		self.ncs_rms = None
		self.nma = None
		self.has_gl = HAS_GL
		self.has_fog = 0
		self.mod = []
		self.nmod = 0
		self.nmol = 0
		self.mol = []
		self.whichmol = None
		self.map = []
		self.sel = []
		self.miniwin = None
		self.renderer = None
		self.rendmod = 3
		self.renwin = None
		self.iren = None
		self.istyle = None
		self.axes = None
		self.camera = None
		self.label = None
		self.islabelon=0
		self.isstarton=0
		self.slabmode=0
		self.fogmode=0
		#~~~~~~Control~~~~~~~~~~~~~~~~~~~~~~~~~~#
		self.slab = 0
		self.slabmax=0
		self.zpos = 0
		self.zposmax = 0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def create_tmp_dir(self):
		if os.path.exists(os.getcwd() + '/veda_tmp'):
			os.system('rm -rf %s'%(os.getcwd() + '/veda_tmp'))
		os.mkdir(os.getcwd() + '/veda_tmp')
		return os.getcwd() + '/veda_tmp'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               Chargement de la fenetre de rendu et de l'interacteur
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def setenv(self,root):
		self.renderer = vtk.vtkOpenGLRenderer()
		#self.renderer.SetBackground(1,1,1)
		self.renwin = vtk.vtkXOpenGLRenderWindow()
		self.renwin.AddRenderer(self.renderer)
		#self.renwin.SetStereoTypeToRedBlue()
		#self.renwin.StereoCapableWindowOn()
		#self.renwin.SetStereoTypeToCrystalEyes()
		#self.renwin.StereoRenderOn() #activer la stereoscopie
		self.renwin.SetStereoTypeToAnaglyph()
		#self.renwin.SetAnaglyphColorMask(4,2) #gauche,droite | r=4, v=2 ,b=1
		self.iren = vtkTkRenderWindowInteractor(root,rw=self.renwin, width=1000, height=800)
		self.istyle = vtk.vtkInteractorStyleSwitch()
		self.iren.SetInteractorStyle(self.istyle)
		self.istyle.SetCurrentStyleToTrackballCamera()
		self.iren.Initialize()
		self.iren.pack(side='top', fill='both', expand=1)
		self.addobservers()
		self.setcam()
		self.setaxes()
		#self.make_fog() #charger le fog par defaut ? /!\ modifier la variable has_fog
		self.iren.Start()
		self.display_ascii_logo()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def embedded_shell_wiz(self):
		gfx=self
		IPython.Shell.hijack_tk()
		ipshell = IPython.Shell.IPShellEmbed(banner= ">>>   [VEDA's Embedded Python Shell]   <<<" )
		ipshell()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def fog_on_off(self):
		if self.has_gl :
			if not self.fogmode:
				self.has_fog=1
				glEnable(GL_FOG)
				glFogf(GL_FOG_MODE, GL_LINEAR)
				glHint(GL_FOG_HINT, GL_NICEST)
				"""
				glFogf(GL_FOG_DENSITY,0.2)
				c=self.camera.GetDistance()
				s=self.slab
				z=self.zpos
				glFogf(GL_FOG_START, c-z-s/2.)
				glFogf(GL_FOG_END, 2000)
				self.renwin.Render()
				"""
				self.update_fog_dist()
				self.update_fog_density()
				self.itf.dist.configure(state=tk.NORMAL)
				self.itf.density.configure(state=tk.NORMAL)
				self.fogmode=1
			else:
				glDisable(GL_FOG)
				self.renwin.Render()
				self.itf.dist.configure(state=tk.DISABLED)
				self.itf.density.configure(state=tk.DISABLED)
				self.fogmode=0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				
	def make_fog(self):
		if self.has_gl :
			self.has_fog=1
			glEnable(GL_FOG)
			glFogf(GL_FOG_MODE, GL_LINEAR)
			glFogf(GL_FOG_DENSITY,0.2)
			glHint(GL_FOG_HINT, GL_NICEST)
			c=self.camera.GetDistance()
			s=self.slab
			z=self.zpos
			glFogf(GL_FOG_START, c-z-s/2.)
			glFogf(GL_FOG_END, c-z+s/2.)
			self.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def disable_fog(self):
			if self.has_fog:
				glDisable(GL_FOG)
				self.renwin.Render()
				self.has_fog=0
			else :
				self.make_fog()				
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def setaxes(self):
		if self.axes==None:
			self.axes = vtk.vtkAxesActor()
			#self.axes.SetShaftTypeToCylinder()
			self.axes.SetTotalLength( 10, 10 , 10 )
			self.axes.SetNormalizedShaftLength( 1, 1, 1 )
			self.axes.SetNormalizedTipLength( 0, 0, 0 )
			#self.axes.SetNormalizedShaftLength( 0.85, 0.85, 0.85 )
			#self.axes.SetNormalizedTipLength( 0.15, 0.15, 0.15 )
			self.axes.AxisLabelsOff()
			axprop = self.axes.GetXAxisTipProperty().SetColor( 0, 0, 1 )
			axprop = self.axes.GetXAxisShaftProperty().SetColor( 0, 0, 1  )
			axprop = self.axes.GetYAxisTipProperty().SetColor( 1, 1, 1 )
			axprop = self.axes.GetYAxisShaftProperty().SetColor( 1, 1, 1 )
			axprop = self.axes.GetZAxisTipProperty().SetColor( 1, 0, 0 )
			axprop = self.axes.GetZAxisShaftProperty().SetColor( 1, 0, 0 )
			self.renderer.AddActor(self.axes)
			self.renwin.Render()
		else :
			if self.axes.GetVisibility():
				self.axes.VisibilityOff()
			else:
				self.axes.VisibilityOn()
			self.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def setcam(self):
		self.camera=vtk.vtkOpenGLCamera()
		self.camera.SetFocalPoint(0, 0, 0)
		self.camera.SetPosition(0, 0, 300)
		self.camera.SetViewUp(0, 0, 0)
		self.camera.SetEyeAngle(5.0)
		self.renderer.SetActiveCamera(self.camera)
		self.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def reset_view(self):
		projmod = self.camera.GetParallelProjection() #parallel=1/perspective=0
		pscale = self.camera.GetParallelScale()
		dist = self.camera.GetDistance()
		self.setcam()
		if projmod:
			self.camera.ParallelProjectionOn()
			self.camera.SetParallelScale(pscale)
		else :
			self.camera.ParallelProjectionOff()
			self.camera.SetDistance(dist)
		self.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def parallel_proj_onoff(self):
		if self.camera.GetParallelProjection():
			 self.camera.ParallelProjectionOff()
		else :
			self.camera.ParallelProjectionOn()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def auto_adjust_clipping_onoff(self):
		if self.istyle.GetAutoAdjustCameraClippingRange():
			self.istyle.AutoAdjustCameraClippingRangeOff()
		else :
			self.istyle.AutoAdjustCameraClippingRangeOn()
		self.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def slab_mode_onoff(self,varslab):
		if self.slabmode==0:
			self.istyle.AutoAdjustCameraClippingRangeOff()
			self.camera.ParallelProjectionOn()
			self.camera.SetParallelScale(self.camera.GetDistance()/4.)
			self.renwin.Render()
			varslab.set(1)
			self.slabmode=1
			self.update_slab()
			self.itf.slab.configure(state=tk.NORMAL)
			self.itf.z.configure(state=tk.NORMAL)
		else :
			self.istyle.AutoAdjustCameraClippingRangeOn()
			self.camera.ParallelProjectionOff()
			self.camera.SetDistance(self.camera.GetParallelScale()*4.)
			self.renderer.ResetCameraClippingRange()
			if self.map!=[]:
				(xmin,xmax,ymin,ymax,zmin,zmax)=self.map[0].box.GetBounds()
				fp=((xmax+xmin)/2.,(ymax+ymin)/2.,(zmax+zmin)/2.)
				self.camera.SetFocalPoint(fp[0],fp[1],fp[2])
			self.renwin.Render()
			varslab.set(0)
			self.slabmode=0
			self.itf.slab.configure(state=tk.DISABLED)
			self.itf.z.configure(state=tk.DISABLED)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def set_open_text_on_off(self):
		if not self.islabelon:
			tm = vtk.vtkTextMapper()
			tp = tm.GetTextProperty()
			tp.SetFontFamilyToArial()
			tp.SetFontSize(20)
			tp.BoldOff()
			tp.ShadowOff()
			tp.SetColor(1, 1, 1)
			tp.SetOpacity(0.8)
			self.label = vtk.vtkActor2D()
			self.label.VisibilityOn()
			self.label.SetMapper(tm)
			v=[200,100]
			self.label.SetPosition(v)
			tm.SetInput("VEDA - Visual Environment for Docking Algorithms - Version (%s)\n \nDevelopped by Gael Goret and Jorge Navaza"%self.version)
			self.renderer.AddActor2D(self.label)
			self.renwin.Render()
			self.islabelon=1
		else :
			self.label.VisibilityOff()
			self.renderer.RemoveActor2D(self.label)
			self.renwin.Render()
			self.islabelon=0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def display_ascii_logo(self):
		ascii_logo=range(11)
		ascii_logo[0]='================================================='
		ascii_logo[1]='  Welcome to  '
		ascii_logo[2]='      __      ________ _____'
		ascii_logo[3]='      \ \    / /  ____|  __ \   /\\'
		ascii_logo[4]='       \ \  / /| |__  | |  | | /  \\'
		ascii_logo[5]='        \ \/ / |  __| | |  | |/ /\ \\'
		ascii_logo[6]='         \  /  | |____| |__| / ____ \\'
		ascii_logo[7]='          \/   |______|_____/_/    \_\\'
		ascii_logo[8]='		             Version : %s'%self.version
		ascii_logo[9]='>> a Visual Environment for Docking Algorithms <<'
		ascii_logo[10]='================================================='
		for i in ascii_logo:
			print i
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def set_open_img_on_off(self):
		if not self.isstarton:
			self.axes.VisibilityOff()
			self.camera.SetFocalPoint(145, 61, 0)
			self.camera.SetPosition(145, 61, 529)
			reader = vtk.vtkTIFFReader()
			reader.SetFileName(self.vedadir + 'doc/veda.tif')
			reader.SetOrientationType(4)
			self.ia = vtk.vtkImageActor()
			self.ia.SetInput(reader.GetOutput())
			self.renderer.AddActor(self.ia)
			self.renwin.Render()
			self.isstarton = 1
		else :
			self.axes.VisibilityOn()
			self.camera.SetFocalPoint(0, 0, 0)
			self.camera.SetPosition(0, 0, 400)
			self.camera.SetViewUp(0, 0, 0)
			self.camera.Elevation(0)
			self.camera.Roll(0)
			self.camera.Azimuth(0)
			self.renderer.RemoveActor(self.ia)
			self.renwin.Render()
			self.isstarton=0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def snapshot(self):
		fn = FD.asksaveasfilename(title="Save snapshot",filetypes=[("JPEG format", "*.jpg "),("TIFF format", "*.tif "),("PNG  format", "*.png ")])
		if fn == '':
			return
		w2i = vtk.vtkWindowToImageFilter()
		w2i.SetInput(self.renwin)
		w2i.Update()
		if fn.endswith('tif') or fn.endswith('TIF') :
			writer = vtk.vtkTIFFWriter()
		elif fn.endswith('jpg') or fn.endswith('JPG') or fn.endswith('jpeg') or fn.endswith('JPEG'):
			writer = vtk.vtkJPEGWriter()
		elif fn.endswith('png') or fn.endswith('PNG'):
			writer = vtk.vtkPNGWriter()
		else:
			fn += '.jpg'
			writer = vtk.vtkJPEGWriter()
		writer.SetFileName(fn)
		writer.SetInputConnection(w2i.GetOutputPort())
		self.renwin.Render()
		writer.Write()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              CONTROL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def addobservers(self):
		""" disable zoom in actor mode but not in camera mode """
		self.istyle.SetCurrentStyleToTrackballActor()
		curst=self.istyle.GetCurrentStyle()
		curst.AddObserver("RightButtonPressEvent",self.rot_in_plane)
		curst.AddObserver("RightButtonReleaseEvent",self.stop_rot_in_plane)
		self.istyle.SetCurrentStyleToTrackballCamera()
		curst=self.istyle.GetCurrentStyle()
		self.iren.bind("<Up>", lambda e, s=self: s.zup(e))
		self.iren.bind("<Down>", lambda e, s=self: s.zdown(e))
		self.iren.bind("<Right>", lambda e, s=self: s.slabup(e))
		self.iren.bind("<Left>", lambda e, s=self: s.slabdown(e))
		#self.iren.bind("<ButtonPress-4>", lambda e, s=self: s.fogup(e))
		#self.iren.bind("<ButtonPress-5>", lambda e, s=self: s.fogdown(e))
		self.iren.bind("<Tab>", lambda e, s=self: s.switch(e))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def switch(self,obj=None,event=None):
		if self.istyle.GetCurrentStyle().GetClassName() =='vtkInteractorStyleTrackballActor':
			self.istyle.SetCurrentStyleToTrackballCamera()
			try:
				cam=tk.PhotoImage(file=self.vedadir + 'doc/cam.gif')
				self.itf.imodebut.configure(image=cam)
				self.itf.imodebut.image=cam
			except : pass
		elif self.istyle.GetCurrentStyle().GetClassName()=='vtkInteractorStyleTrackballCamera':
			self.istyle.SetCurrentStyleToTrackballActor()
			try :
				mol=tk.PhotoImage(file=self.vedadir + 'doc/mol.gif')
				self.itf.imodebut.configure(image=mol)
				self.itf.imodebut.image=mol
			except : pass
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def rot_in_plane(self,obj,event):
		self.iren.SetControlKey(1)
		obj.OnLeftButtonDown()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def stop_rot_in_plane(self,obj,event):
		self.iren.SetControlKey(0)
		obj.OnLeftButtonUp()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def fogup(self,event):
		if self.zpos < self.zmax-3 :
			self.zpos += 3
			self.update_fog()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def fogdown(self,event):
		if self.zpos > 3 -self.zmax :
			self.zpos -= 3
			self.update_fog()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def slabup(self,event):
		""" increase slab size"""
		if self.slab < self.slabmax-2 :
			self.itf.slabvar.set(self.itf.slabvar.get()+2)
			self.slab += 2
			self.update_slab()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def slabdown(self,event):
		"""decrease slab size """
		if self.slab > 2 :
			self.itf.slabvar.set(self.itf.slabvar.get()-2)
			self.slab -= 2
			self.update_slab()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def slabchange(self,event):
		self.slab = self.itf.slabvar.get()
		self.update_slab()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def zup(self,event):
		"""decrease z-position of slab down """
		if self.zpos > 2-self.zposmax :
			self.zpos -= 2
			self.itf.zvar.set(self.itf.zvar.get()-2)
			self.update_slab()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def zdown(self,event):
		""" increase z-position of slab """
		if self.zpos < self.zposmax-2 :
			self.itf.zvar.set(self.itf.zvar.get()+2)
			self.zpos += 2
			self.update_slab()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def zchange(self,event):
		self.zpos = self.itf.zvar.get()
		self.update_slab()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def update_slab(self):
		""" changes slab of rendering"""
		c = self.camera.GetDistance()
		s = self.slab
		z = self.zpos
		self.camera.SetClippingRange(c-z-s/2.,c-z+s/2.)
		self.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def update_fog_dist(self,obj=None,event=None):
		""" changes fog of rendering"""
		c=self.camera.GetDistance()
		fd=self.itf.fogdivar.get()
		if self.has_gl:
			glFogf(GL_FOG_START,c-fd)
			glFogf(GL_FOG_END, 5*c)
			self.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def update_fog_density(self,obj=None,event=None):
		""" changes fog of rendering"""
		dens=self.itf.fogdevar.get()/100.
		if self.has_gl:
			glFogf(GL_FOG_DENSITY,dens)
			self.renwin.Render()
