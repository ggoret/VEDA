# -*- coding: utf-8 -*-
# NMA Module of VEDA
#
# Copyright:
# 2009-01-01 Gael Goret  gael.goret@ibs.fr
#
# Last modified:
# 2010-12-10 Gael Goret  gael.goret@ibs.fr
# deprecate (directory changed)
try :
	import Tkinter as tk
except :
	import tkinter as tk
import vtk
from vtk.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor
import gfx ,Map,mod, itf,sym,uro
from os import system,chdir,path,mkdir
from numpy import *
import time
import tkMessageBox as MB
import tkFileDialog as FD
class normal_modes():
	def __init__(self,gfx):
		self.gfx=gfx
		self.mod=[]
		self.cutoff=15. #CUTOFF
		self.force_const=10. #FORCE CONST
		self.mode=None
		self.rms=None
		self.nbpas=None
		self.viewer=None
		self.currentmol=None
		self.nbstep=None
		self.oldon=None
		self.loop=None
		self.loopobs=None
		gfx.nma = self	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def create_mod_dir(self,gfx,mod):
		if not path.exists(gfx.tmpdir + '/%s'%mod.un):
			mkdir(gfx.tmpdir + '/%s'%mod.un)	
		mod.nmfn =gfx.tmpdir + '/%s/'%mod.un + extract_file_from_path(mod.rfn)[:-6] + '_CA.pdb'
		system("grep ' CA ' %s | grep '^ATOM' > %s"%(mod.rfn,mod.nmfn))
		try:
		   	dat = open(gfx.tmpdir + '/%s/'%mod.un + '/pdbmat.dat','w')
		except:
		   	print 'Problem building pdbmat.dat'
			return
		dat.write('Coordinate FILENAME        = %s\n'%extract_file_from_path(mod.nmfn))
		dat.write('INTERACtion DISTance CUTOF =      %s\n'%self.cutoff)
		dat.write('INTERACtion FORCE CONStant =      %s\n'%self.force_const)
		dat.write('Origin of MASS values      =       CONS ! CONstant, or from COOrdinate file.\n')
		dat.write('Output PRINTing level      =          1 ! =1: more detailled. =2: debug level.\n')
		dat.write('Bond DEFINITION            =       NONE ! NONe, ALL, or between CONsecutive atoms.\n')
		dat.write('Maximum bond LENGTH        =      0.000\n')
		dat.write('LevelSHIFT                 =    1.0E-09 ! Non-zero value often required (numerical reasons).\n')
		dat.write('Matrix FORMAT              =     FREE   ! Free, or Binary, matrix saved.\n')
		dat.close()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				
	def build_hessian_mat(self,gfx,mod):
		chdir(gfx.tmpdir + '/%s/'%mod.un)
		gfx.itf.status.set('computation of the hessian in progress...')
		system(gfx.vedabin + '/pdbmat.exe')
		gfx.itf.status.clear()
		chdir(gfx.workdir)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def diag_hessian_mat(self,gfx,mod):
		chdir(gfx.tmpdir + '/%s/'%mod.un)	
		gfx.itf.status.set('Diagonalisation in progress ... can take few (or so much) minutes ')
		system(gfx.vedabin + '/diagstd.exe &')
		gfx.itf.status.clear()
		chdir(gfx.workdir)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def shift(self,gfx,mod,mode,rms,nbpas):
		chdir(gfx.tmpdir + '/%s/'%mod.un)
		if path.exists('out1'):
			system('rm -f out*')
		gfx.itf.status.set('Amplitude sampling in progress ...')
		system(gfx.vedabin + '/nmashift.sh %s %s %s %s'%(extract_file_from_path(mod.nmfn),mode,rms,nbpas))
		gfx.itf.status.clear()
		chdir(gfx.workdir)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def display_conformers_manual(self,obj=None,event=None):
		if self.currentmol!=None:
			arg = self.gfx.itf.varconformerSitu.get()+self.nbstep
			if self.oldon==None:
				self.currentmol.nmcl[0].VisibilityOff()
				self.currentmol.nmcl[arg].VisibilityOn()
				self.gfx.renwin.Render()
				self.oldon=arg
			else :
				self.currentmol.nmcl[self.oldon].VisibilityOff()
				self.currentmol.nmcl[arg].VisibilityOn()
				self.oldon=arg
				self.gfx.renwin.Render()
			self.gfx.renderer.ResetCameraClippingRange()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def display_in_situ(self,gfx,mol,nbstep):
		gfx.itf.status.set('Loading all conformers, progression : 0.0 percents ')
		spectra=[]
		cpt=0
		self.nbstep=nbstep
		for outnb in range(-nbstep,nbstep+1):
			fn= gfx.tmpdir + '/%s/'%mol.mod.un + 'out%s'%(outnb)
			if path.exists(fn):
				cpt+=1
				pdb = open(fn,'r')
				nat=0
				vx=[]
				vy=[]
				vz=[]
				for l in pdb:
					if l[0:4] == 'ATOM':
						nat+=1
						x=float(l[30:38])
						y=float(l[38:46])
						z=float(l[46:54])
						vx += [x]
						vy += [y]
						vz += [z]
				pdb.close()
				camax = 50
				reader=vtk.vtkAppendPolyData()
				for i in range (nat-1) :
				    	#if ( (vx[i]-vx[i+1])**2+(vy[i]-vy[i+1])**2+(vz[i]-vz[i+1])**2 ) < camax :#avoid long lines
			    		s = vtk.vtkLineSource()
					s.SetPoint1([vx[i],vy[i],vz[i]])
					s.SetPoint2([vx[i+1],vy[i+1],vz[i+1]])
					reader.AddInput(s.GetOutput())
				mapper = vtk.vtkPolyDataMapper()
				mapper.SetInputConnection(reader.GetOutputPort())
				mapper.UseLookupTableScalarRangeOff()
				mapper.SetScalarVisibility(1)
				mapper.SetScalarModeToDefault()
				acteur = vtk.vtkOpenGLActor()
				acteur.SetMapper(mapper)
				if mol.acteur.GetClassName()=='vtkAssembly':
					acteur.GetProperty().SetColor(mol.acteur.GetParts().GetItemAsObject(0).GetProperty().GetColor())
					acteur.GetProperty().SetLineWidth(mol.acteur.GetParts().GetItemAsObject(0).GetProperty().GetLineWidth())
				else : 
					acteur.GetProperty().SetColor(mol.acteur.GetProperty().GetColor())
					acteur.GetProperty().SetLineWidth(mol.acteur.GetProperty().GetLineWidth())
				assigne_rotra(mol,acteur)
				spectra+=[acteur]
				gfx.renderer.AddActor(acteur)
				acteur.VisibilityOff()
				gfx.itf.status.set('Loading all conformers, progression : %3.1f percents'%((100.*cpt)/(nbstep*2.+1)))
				gfx.renwin.Render()
		gfx.renderer.ResetCameraClippingRange()
		gfx.renwin.Render()
		mol.nmcl = spectra[:]
		self.currentmol=mol
		self.oldon=None
		self.currentmol.nmcl[0].VisibilityOn()
		self.currentmol.acteur.VisibilityOff()
		if self.currentmol.lsm != []:
			self.currentmol.lsm[0].VisibilityOff()
		gfx.renwin.Render()
		gfx.renderer.ResetCameraClippingRange()
		self.loop = Callbackinsitu(gfx)
		self.loopobs = gfx.iren.AddObserver('TimerEvent', self.loop.loop)
		gfx.iren.CreateRepeatingTimer(100)
		gfx.iren.TimerEventResetsTimerOn()
		gfx.iren.Start()
		gfx.itf.status.clear()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def stop_insitu_loop(self,gfx):
		if self.loop!=None:
			gfx.iren.TimerEventResetsTimerOff()
			gfx.iren.RemoveObserver(self.loopobs)
			self.currentmol.acteur.VisibilityOff()
			self.loop=None
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def reset_normal_state(self,gfx):
		if self.currentmol!=None and self.oldon!=None:
			self.currentmol.nmcl[self.oldon].VisibilityOff()
			self.currentmol.acteur.VisibilityOn()
		gfx.renwin.Render()
		gfx.renderer.ResetCameraClippingRange()	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def save_conf_as(self,gfx):
		conf = self.currentmol.nmcl[self.gfx.itf.varconformerSitu.get()+self.nbstep]
		fn = FD.asksaveasfilename(parent=self.gfx.itf.root,filetypes=[('PDB files','*.pdb')] ,title="Save conformer %d coordinates as..."%(self.gfx.itf.varconformerSitu.get()+self.nbstep))
		if len(fn ) > 0:
			if not fn.endswith('.pdb'):
				fn = fn +'.pdb'
    			mat = conf.GetMatrix() #recuperation de la mat de pv 4x4
    			R = array(sym.ExtractRotMat(mat)) #extraction de la matrice de rotation depuit la vtk 4x4
			t = array([mat.GetElement(0,3),mat.GetElement(1,3),mat.GetElement(2,3)]) #extract tx,ty,tz
			pdbin=open(gfx.tmpdir + '/%s/'%self.currentmol.mod.un + 'out%s'%(self.gfx.itf.varconformerSitu.get()),"r")
			pdbout=open(fn,"w")
			for ligne in pdbin:
				if (ligne[0:6] == "ATOM  " or ligne[0:6] == "HETATM"):
					x = float(ligne[30:38])
					y = float(ligne[38:46])
					z = float(ligne[46:54])
					r = array([x,y,z])
					newr = dot(R,r) + t
					ligne=ligne[:30]+"%8.3f%8.3f%8.3f"%tuple(newr)+ligne[54:]
				pdbout.write(ligne)
			pdbin.close()
			pdbout.close()
			return
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def norma(self,gfx,lowstep,highstep):
		try:
			lowstep=int(lowstep)
			highstep=int(highstep)
		except:
			MB.showwarning('Warning','Input Problem for Interval variables')
			return
		self.build_nm_fitin_file()
		if gfx.mol!=[]:
			if gfx.fit!=None and uro.fit_is_possible(gfx):
				if uro.sym_is_ok(gfx):
					system('e/norma %s %d %d &'%(self.currentmol.mod.un,lowstep,highstep))
					gfx.itf.showinfo('Info','URO Refinement is running in background for the different conformers')
				else :
					MB.showwarning('Info','There is no sym-mate molecule within the volume. Update Sym-Mates')
			else :
				MB.showwarning('Info','Configure Fitting|Setup')
		else :
			MB.showwarning('Info','Go to Assemblage|Molecules and load molecules. Then configure Fitting|Setup')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def build_nm_fitin_file(self):
		no = open('N0','w')
		for mol in self.gfx.mol:
			mat=mol.acteur.GetMatrix()
			nbsymop = len(mol.lnbsm)
			strsymlist = str(mol.lnbsm)[1:-1].replace(',','')
			(a,b,g) = sym.R2Eul(sym.ExtractRotMat(mat))
			(x,y,z) = [mat.GetElement(0,3),mat.GetElement(1,3),mat.GetElement(2,3)]
			if mol.id == self.currentmol.id:
				cur = ' #%2d %.3f %.3f %.3f %.3f %.3f %.3f %d %s\n'%(33,a,b,g,x,y,z,nbsymop,strsymlist)
			else :
				no.write(' #%2d %.3f %.3f %.3f %.3f %.3f %.3f %d %s\n'%(mol.mod.id,a,b,g,x,y,z,nbsymop,strsymlist))
				
		no.write(cur)
		no.close()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       		Fonctions Hors Class
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_nm(gfx,modlist3,nmmodlist,nmmollist):
	if gfx.mod!=None:
		for mod in gfx.mod:
			if modlist3.get(tk.ACTIVE) == mod.un:
				if mod.type=='mol':
					if gfx.nma==None:
						nm=normal_modes(gfx)
					if mod in gfx.nma.mod:
						MB.showwarning('Info','The normal modes for this model already exists')
						return
					gfx.nma.create_mod_dir(gfx,mod)
					gfx.nma.build_hessian_mat(gfx,mod)
					gfx.nma.diag_hessian_mat(gfx,mod)
					gfx.nma.mod += [mod]
					build_lists(gfx,nmmodlist,nmmollist)
				else:
					MB.showwarning('Info','Normal modes analysis for map models not yet implemented')
					return
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_nm(gfx,modlist3,nmmodlist,nmmollist):
	if gfx.mod!=None:
		for mod in gfx.mod:
			if modlist3.get(tk.ACTIVE) == mod.un:
				if mod.type=='mol':
					if gfx.nma==None:
						nm=normal_modes(gfx)
					if mod in gfx.nma.mod:
						MB.showwarning('Info','The normal modes for this model already exists')
						return
					fn = FD.askopenfilename(title="Open NM Matrix file",filetypes=[("NM Matrix", "*.eigenfacs"),("All", "*.*")])
					if len(fn)>0:
						if not path.exists(gfx.tmpdir + '/%s'%mod.un):
							system('mkdir %s'%(gfx.tmpdir + '/%s'%mod.un))
						system('cp -f %s %s'%(fn,(gfx.tmpdir + '/%s/matrix.eigenfacs'%mod.un)))
						mod.nmfn =gfx.tmpdir + '/%s/'%mod.un + extract_file_from_path(mod.rfn)[:-6] + '_CA.pdb'
						system("grep ' CA ' %s | grep '^ATOM' > %s"%(mod.rfn,mod.nmfn))
						gfx.nma.mod += [mod]
						build_lists(gfx,nmmodlist,nmmollist)
				else:
					MB.showwarning('Info','Normal modes analysis for map models not yet implemented')
					return
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def save_nm(gfx,nmmodlist):
	if gfx.nma!=None:
		for mod in gfx.nma.mod:
			if nmmodlist.get(tk.ACTIVE) == mod.un:
				fn = FD.asksaveasfilename(parent=gfx.itf.root,filetypes=[('NM Matrix','*.eigenfacs')] ,title="Save normal modes matrix of model %s as..."%mod.un)
				if len(fn) > 0:
					if not fn.endswith('.eigenfacs'):
						fn = fn +'.eigenfacs'
					system('cp -f %s %s'%((gfx.tmpdir + '/%s/matrix.eigenfacs'%mod.un),fn))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def remove_mod(gfx,nmmodlist,nmmollist):
	if gfx.nma!=None:
		for mod in gfx.nma.mod:
			if nmmodlist.get(tk.ACTIVE) == mod.un:
				if path.exists(gfx.tmpdir + '/%s'%mod.un):
					system('rm -rf %s'%(gfx.tmpdir + '/%s'%mod.un))	
				for i in range(len(gfx.nma.mod)):
					if gfx.nma.mod[i].un == mod.un:
						del gfx.nma.mod[i]
				build_lists(gfx,nmmodlist,nmmollist)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def remove_nmmod(gfx,model):
	if gfx.nma!=None:
		for mod in gfx.nma.mod:
			if model.un == mod.un:
				if path.exists(gfx.tmpdir + '/%s'%mod.un):
					system('rm -rf %s'%(gfx.tmpdir + '/%s'%mod.un))	
				for i in range(len(gfx.nma.mod)):
					if gfx.nma.mod[i].un == mod.un:
						del gfx.nma.mod[i]
				if gfx.itf.nmwizopen:
					build_lists(gfx,gfx.itf.nmmodlist,gfx.itf.nmmollist)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def display_nm(gfx,mod,mode,rms,nbpas):
	mode = int(mode)
	rms = float(rms)
	nbpas = int(nbpas)
	gfx.nma.shift(gfx,mod,mode,rms,nbpas)
	gfx.nma.viewer.load_spectra(gfx,mod,nbpas)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def replace_in_situ(gfx,mol,mode,rms,nbpas):
	if not matrix_is_ready(gfx,mol.mod):
		MB.showwarning('Info','Normal Modes are still computing ... please wait')
		return
	mode = int(mode)
	if mode<7 or mode >26:
		MB.showwarning('Error','Incorrect mode number, must be into [7..26] ')
		return
	rms = float(rms)
	nbpas = int(nbpas)
	if gfx.nma.currentmol != mol or gfx.nma.mode != mode or gfx.nma.rms != rms or gfx.nma.nbpas != nbpas:
		gfx.nma.mode=mode
		gfx.nma.rms=rms
		gfx.nma.nbpas=nbpas
		if gfx.nma.loop != None:
			gfx.nma.stop_insitu_loop(gfx)
			gfx.nma.reset_normal_state(gfx)
		#test here if size of matrix != 0
		gfx.nma.shift(gfx,mol.mod,mode,rms,nbpas)
		gfx.nma.display_in_situ(gfx,mol,nbpas)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def matrix_is_ready(gfx,mod):
	mfile  = gfx.tmpdir + '/%s/matrix.eigenfacs'%mod.un 
	if path.getsize(mfile)==0:
		return 0
	else : 
		return 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Callbackinsitu():
	def __init__(self,gfx):
		self.gfx=gfx
		self.sens=1
		self.currentconf=None
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
	def loop(self,obj,event):
		nbc = len(self.gfx.nma.currentmol.nmcl)
		if self.currentconf!=None:
			i=self.currentconf
		else : 
			i = 2
		if i==1 or i==nbc:
			self.sens = self.sens * -1
		i = i + self.sens
		self.currentconf=i
		conf=i-1
		if self.gfx.nma.oldon==None:
			self.gfx.nma.currentmol.nmcl[0].VisibilityOff()
			self.gfx.nma.currentmol.nmcl[conf].VisibilityOn()
			self.gfx.renwin.Render()
			self.gfx.nma.oldon=conf
		else :
			self.gfx.nma.currentmol.nmcl[self.gfx.nma.oldon].VisibilityOff()
			self.gfx.nma.currentmol.nmcl[conf].VisibilityOn()
			self.gfx.nma.oldon=conf
			self.gfx.renwin.Render()
		self.gfx.renderer.ResetCameraClippingRange()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       		MISC functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			
def build_lists(gfx,nmmodlist,nmmollist):
	if gfx.itf.nmwizopen:
		nmmodlist.delete(0, tk.END)
		nmmollist.delete(0, tk.END)
		if gfx.nma!=None:
			for mod in gfx.nma.mod:
				nmmodlist.insert(tk.END,mod.un)
			for mol in gfx.mol:
				for mod in gfx.nma.mod:
					if mol.mod.un == mod.un:
						nmmollist.insert(tk.END,mol.un)
						nmmollist.itemconfigure(tk.END,background=vtk2tkhex_color(mol.col))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def vtk2tkhex_color(vtk_col):
	(r,v,b)=(vtk_col[0]*255.0, vtk_col[1]*255.0, vtk_col[2]*255.0)
	tk_rgb = "#%02x%02x%02x" % (r, v, b)
  	return tk_rgb
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def extract_file_from_path(fn):
	for i in range(len(fn)):
		if fn[-(i+1)]=='/':
			return fn[-i:]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def assigne_rotra(mol,confact):
	confact.SetPosition(mol.acteur.GetPosition())
	confact.SetOrientation(mol.acteur.GetOrientation())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       		Normal Modes Viewer
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class nmviewer():
	def __init__(self,gfx,frame):
		self.renderer = vtk.vtkOpenGLRenderer()
		self.renwin = vtk.vtkXOpenGLRenderWindow()
		self.renwin.AddRenderer(self.renderer)
		self.renwin.SetStereoTypeToAnaglyph()
		self.iren = vtkTkRenderWindowInteractor(frame,rw=self.renwin, width=500, height=500 )
		self.istyle = vtk.vtkInteractorStyleSwitch()
		self.iren.SetInteractorStyle(self.istyle)
		self.istyle.SetCurrentStyleToTrackballCamera()
		self.iren.Initialize()
		self.iren.pack(side='bottom', fill='both', expand=0)
		self.iren.Start()
		self.camera=vtk.vtkCamera()
		self.camera.SetFocalPoint(0, 0, 0)
		self.camera.SetPosition(0, 0, 250)
		self.camera.SetViewUp(0, 0, 0)
		self.camera.SetEyeAngle(5.0)
		self.tm=None
		self.renderer.SetActiveCamera(self.camera)
		self.renwin.Render()
		self.currentmod= None
		self.nbstep=None
		self.oldon=None
		self.tloop=None
		self.tloopobs=None
		self.gfx=gfx
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def display_progression(self):
		self.tm = vtk.vtkTextMapper()
		tp = self.tm.GetTextProperty()
		#tp.SetFontFamilyToArial()
		tp.SetFontSize(25)
		#tp.BoldOff()
		tp.ShadowOff()
		tp.SetColor(1, 1, 1)
		tp.SetOpacity(1)
		label = vtk.vtkActor2D()
		label.VisibilityOn()
		label.SetMapper(self.tm)
		v=[20,20]
		label.SetPosition(v)
		self.tm.SetInput("Loading, please wait ... ")
		self.renderer.AddActor2D(label)
		self.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def load_spectra(self,gfx,mod,nbstep):
		self.renderer.RemoveAllViewProps()
		gfx.itf.status.set('Loading all conformers, please wait ...')
		spectra=[]
		cpt=0
		self.nbstep=nbstep
		self.display_progression()
		for outnb in range(-nbstep,nbstep+1):
			fn= gfx.tmpdir + '/%s/'%mod.un + 'out%s'%(outnb)
			if path.exists(fn):
				cpt+=1
				pdb = open(fn,'r')
				nat=0
				vx=[]
				vy=[]
				vz=[]
				for l in pdb:
					if l[0:4] == 'ATOM':
						nat+=1
						x=float(l[30:38])
						y=float(l[38:46])
						z=float(l[46:54])
						vx += [x]
						vy += [y]
						vz += [z]
				pdb.close()
				camax = 50
				reader=vtk.vtkAppendPolyData()
				for i in range (nat-1) :
				    	#if ( (vx[i]-vx[i+1])**2+(vy[i]-vy[i+1])**2+(vz[i]-vz[i+1])**2 ) < camax :#avoid long lines
			    		s = vtk.vtkLineSource()
					s.SetPoint1([vx[i],vy[i],vz[i]])
					s.SetPoint2([vx[i+1],vy[i+1],vz[i+1]])
					reader.AddInput(s.GetOutput())
				mapper = vtk.vtkPolyDataMapper()
				mapper.SetInputConnection(reader.GetOutputPort())
				mapper.UseLookupTableScalarRangeOff()
				mapper.SetScalarVisibility(1)
				mapper.SetScalarModeToDefault()
				acteur = vtk.vtkOpenGLActor()
				acteur.SetMapper(mapper)
				acteur.GetProperty().SetColor(0,1,0)
				acteur.GetProperty().SetLineWidth(3)
				acteur.SetMapper(mapper)
				spectra+=[acteur]
				self.renderer.AddActor(acteur)
				acteur.VisibilityOff()
				self.tm.SetInput("%3.1f"%((100.*cpt)/(nbstep*2.+1))+' %')
				self.renwin.Render()
		self.tm.SetInput("")
		self.renderer.ResetCameraClippingRange()
		self.renwin.Render()
		mod.nmcl = spectra[:]
		self.currentmod=mod
		self.currentmod.nmcl[self.nbstep].VisibilityOn()#on allume la structure non déformée
		self.oldon=None
		self.renwin.Render()
		self.renderer.ResetCameraClippingRange()
		gfx.itf.confscroll.configure(from_=-self.nbstep,to=self.nbstep)
		gfx.itf.varconformer.set(0)
		gfx.itf.status.clear()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def display_conformers_manual(self,obj=None,event=None):
		if self.currentmod!=None:
			arg = self.gfx.itf.varconformer.get()+self.nbstep
			if self.oldon==None:
				self.currentmod.nmcl[self.nbstep].VisibilityOff()
				self.currentmod.nmcl[arg].VisibilityOn()
				self.renwin.Render()
				self.oldon=arg
			else :
				self.currentmod.nmcl[self.oldon].VisibilityOff()
				self.currentmod.nmcl[arg].VisibilityOn()
				self.oldon=arg
				self.renwin.Render()
			self.renderer.ResetCameraClippingRange()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def loop_thread(self):
		if self.currentmod != None:
			if self.gfx.itf.varloop.get()==1:
	   			self.iren.ReInitialize()
				self.tloop = Callbackviewer(self,self.gfx)
	   			self.tloopobs = self.iren.AddObserver('TimerEvent', self.tloop.loop)
	   			self.tlooptimerid=self.iren.CreateRepeatingTimer(100)
	   			self.iren.TimerEventResetsTimerOn()
				self.iren.Start()
			else :
				self.iren.TimerEventResetsTimerOff()
				self.iren.RemoveObserver(self.tloopobs)
		else : 
			MB.showwarning('Info','Click on Display button first')
			self.gfx.itf.varloop.set(0)
			return
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Callbackviewer():
	def __init__(self,nmviewer,gfx):
		self.viewer=nmviewer
		self.gfx=gfx
		self.sens=1
		self.currentconf=None
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	def loop(self,obj,event):
		nbc = len(self.viewer.currentmod.nmcl)
		if self.currentconf == None:
			i=self.gfx.itf.varconformer.get()+self.viewer.nbstep+1
			if self.gfx.itf.varconformer.get() < 0:
				self.sens=-1
			else :
				self.sens=1
		else :
			i= self.currentconf
		if i==1 or i==nbc:
			self.sens = self.sens * -1
		i = i + self.sens
		self.currentconf=i
		conf=i-1
		if self.viewer.currentmod!=None:
			self.gfx.itf.varconformer.set(conf-self.viewer.nbstep)
		if self.viewer.oldon==None:
			self.viewer.currentmod.nmcl[self.viewer.nbstep].VisibilityOff()
			self.viewer.currentmod.nmcl[conf].VisibilityOn()
			self.viewer.renwin.Render()
			self.viewer.oldon=conf
		else :
			self.viewer.currentmod.nmcl[self.viewer.oldon].VisibilityOff()
			self.viewer.currentmod.nmcl[conf].VisibilityOn()
			self.viewer.oldon=conf
			self.viewer.renwin.Render()
		self.viewer.renderer.ResetCameraClippingRange()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################################END##########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
