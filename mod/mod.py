# -*- coding: utf-8 -*-
# Models and Molecules Module of VEDA
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
import gfx,Map,itf,sym,uro,Nma
from os import system,chdir,path
import tkColorChooser
import tkMessageBox as MB
import tkFileDialog as FD
from numpy import *
from vtk.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor
#===============================================================================
#                              Definition du modele + wizard function
#===============================================================================
class Mod():
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def __init__(self):
		self.type=None# 'map' or 'mol'
		self.dspmodtype=None# 'AA' or 'B' or 'CA'
		self.id=None #unique identifiant
		self.fn=None #filename
		self.dfn=None#display file name
		self.rfn=None#reference file name
		self.ofn=None#origin file name
		self.nmfn=None#nm file name
		self.un=None #unique name
		self.nmol=0
		self.rottra=None#rotation and translation
		self.coldist=None#collision distance
		self.sphdist=None
		self.rep=None# representation type
		self.sigavg=None #sigma and avg value for map
		self.isov=None
		self.outputfile=None
		self.nmcl=[]#Normal modes conformeres list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def assign_mod(gfx,nfv,cpltnfv,nameentry):
	cpltnfv.set(itf.browsemol())
	nfv.set(extract_file_from_path(cpltnfv.get()))
	nameentry.delete(0,tk.END)
	nameentry.insert(tk.END,'M%d'%(gfx.nmod+1))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def add_mod(gfx,modlist,nameentry,nfv,cpltnfv,varrepsv,isov):
	modname=nameentry.get()
	varrep=varrepsv.get() #stringvar
	if cpltnfv.get() == u'':
		MB.showwarning('Info','Select a file')
		return
	if modname == '' :
		MB.showwarning('Info','Name entry empty, fill in the blank')
		return
	for model in gfx.mod:
		if model.un == modname :
			MB.showwarning('Info','Name already defined, chose another one')
			return
	if varrep=='':
		MB.showwarning('Info','Select a model type')
		return
	gfx.itf.root.configure(cursor='watch')
	gfx.itf.status.set('Model loading ... please wait')
	modo=Mod()
	modo.id=set_mod_id(gfx)
	modo.ofn = fn = cpltnfv.get()
	if varrep == 'Backbone' or varrep =='All Atoms' or varrep =='C alpha':#PDB CASE
		rfn =gfx.tmpdir + '/' + extract_file_from_path(fn)[:-4] + '_R.pdb'
		bfn =gfx.tmpdir + '/' + extract_file_from_path(fn)[:-4] + '_B.pdb'
		cafn =gfx.tmpdir + '/' + extract_file_from_path(fn)[:-4] + '_CA.pdb'
		if path.exists(rfn):
			i=1
			while path.exists(gfx.tmpdir + '/' + extract_file_from_path(fn)[:-4] + '_R_%d.pdb'%i):
				i+=1
			rfn = gfx.tmpdir + '/' + extract_file_from_path(fn)[:-4] + '_R_%d.pdb'%i
		if path.exists(bfn):
			i=1
			while path.exists(gfx.tmpdir + '/' + extract_file_from_path(fn)[:-4] + '_B_%d.pdb'%i):
				i+=1
			bfn = gfx.tmpdir + '/' + extract_file_from_path(fn)[:-4] + '_B_%d.pdb'%i
		if path.exists(cafn):
			i=1
			while path.exists(gfx.tmpdir + '/' + extract_file_from_path(fn)[:-4] + '_CA_%d.pdb'%i):
				i+=1
			cafn = gfx.tmpdir + '/' + extract_file_from_path(fn)[:-4] + '_CA_%d.pdb'%i
		chdir(gfx.tmpdir)
		system(gfx.vedabin + '/refer.exe <<ENDOF\n %s  %s \nENDOF'%(fn,rfn))
		if isca(modo.ofn) and (varrep == 'Backbone' or varrep =='All Atoms'):
			varrep='C alpha'
			varrepsv.set('C alpha')
			MB.showwarning('Warning','Your model contains only C-alpha')
		if varrep=='C alpha':
			if not is_type_consistent(gfx,'CA'):
				gfx.itf.status.clear()
				gfx.itf.root.configure(cursor='arrow')
				return
			modo.dfn=cafn
			modo.dspmodtype = 'CA'
			system("grep ' CA ' %s | grep '^ATOM' > %s"%(rfn,modo.dfn))
		elif varrep=='Backbone':
			if not is_type_consistent(gfx,'B'):
				gfx.itf.status.clear()
				gfx.itf.root.configure(cursor='arrow')
				return
			modo.dfn=bfn
			modo.dspmodtype = 'B'
			extract_backbone(rfn,modo.dfn)
		else :
			if not is_type_consistent(gfx,'AA'):
				gfx.itf.status.clear()
				gfx.itf.root.configure(cursor='arrow')
				return
			modo.dspmodtype = 'AA'
			modo.dfn=rfn
		chdir(gfx.workdir)
		modo.rottra = get_rottra(gfx,rfn)
		modo.fn = modo.rfn = rfn
		modo.type='mol'
	elif varrep == 'Surface': #MAP CASE
		if  fn.endswith('.ezd'):
			chdir(gfx.tmpdir)
			mapfileout = extract_file_from_path(fn)[:-4]+'_R.vtk'
			if path.exists(gfx.tmpdir + '/' + mapfileout):
				i=1
				while path.exists(gfx.tmpdir + '/' + mapfileout[:-6] + '_R_%d.vtk'%i):
					i+=1
				mapfileout = mapfileout[:-6] + '_R_%d.vtk'%i
			e2v_out='info_map_%s'%modname
			system(gfx.vedabin+'/e2v.exe >> %s <<ENDOF\n%s  \n%f  \n%s  \nENDOF'%(e2v_out,fn,1,mapfileout)) #1=scale
			modo.dfn = modo.rfn = gfx.tmpdir + '/' + mapfileout
			modo.sigavg = Map.map_sigma_avg(e2v_out)
			(ox,oy,oz)= import_origin(e2v_out)
			system("sed -i -e /^ORIGIN/s/.*/'ORIGIN %f %f %f'/ %s"%(ox,oy,oz,modo.dfn))
			modo.isov = float(isov)
			chdir(gfx.workdir)
		modo.fn=fn
		modo.type='map'
	nfv.set(extract_file_from_path(modo.rfn)) #Reference file name
	modo.rep=varrep
	modo.un=modname
	gfx.mod+=[modo]
	refresh_modlist(gfx,modlist)
	if gfx.itf.molwizopen:
			refresh_modlist(gfx,gfx.itf.modlist2)
	if gfx.itf.nmwizopen:
			refresh_modlist(gfx,gfx.itf.modlist3)
	gfx.miniwin.display_mod(modo)
	gfx.itf.status.clear()
	gfx.itf.root.configure(cursor='arrow')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def is_type_consistent(gfx,currenttype):
	for mod in gfx.mod:
		if currenttype != mod.dspmodtype:
			return MB.askyesno('Warning', "The model type is different to the previous ones. this will cause problem during fitting, would you like to continue ?")
	return 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def extract_backbone(rfn,dfn):
	f=open(rfn,'r')
	fp=open(dfn,'w')
	for l in f:
		if l[13:15] in ['CA','N ','C ']:
			fp.write(l)
	f.close()
	fp.close()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def import_origin(info_map) :
	try:
		f=open(info_map,"r")
	except IOError:
		MB.showwarning('Warning','Conversion from ezd to vtk failed, file info_map_* does not exist')
		return
	for l in f:
		if 'ORIGIN' in l:
			ox,oy,oz = float(l.split()[-3]),float(l.split()[-2]),float(l.split()[-1])
	f.close()
	return (ox,oy,oz)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def extract_file_from_path(fn):
	for i in range(len(fn)):
		if fn[-(i+1)]=='/':
			return fn[-i:]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_rottra(gfx,fn):
	pdb = open(fn,'r')
	rottra=[]
        for l in pdb:
		if l[0:6] == 'ROTTRA':
			lst = l.split()
			rottra = [ float(lst[1]), float(lst[2]), float(lst[3]), float(lst[4]), float(lst[5]), float(lst[6]) ]
	pdb.close()
	return rottra
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def remove_mod(gfx,modlist,modnameentry,nfv):
	if modlist.curselection()!=():
		gfx.miniwin.renderer.RemoveAllViewProps()
		gfx.miniwin.renwin.Render()
	if gfx.mod != []:
		for selmod in modlist.curselection():
			for i in range(len(gfx.mod)):
				if modlist.get(selmod)==gfx.mod[i].un:
					nmolformod = 0
					for mol in gfx.mol:
						if mol.mod.id==gfx.mod[i].id:
							nmolformod+=1
					while nmolformod!=0:
						for j in range(len(gfx.mol)):
							if gfx.mod[i].id==gfx.mol[j].mod.id:
								for sm in gfx.mol[j].lsm:
									gfx.renderer.RemoveActor(sm)
									del sm
								gfx.renderer.RemoveActor(gfx.mol[j].acteur)
								if gfx.mol[j].lsphere!=None:
									gfx.renderer.RemoveActor(gfx.mol[j].lsphere)
								unkeys = gfx.itf.mollist.get(0,tk.END)
								for k in range(len(unkeys)):
									if unkeys[k] == gfx.mol[j].un:
										gfx.itf.mollist.delete(k)
										break
								del gfx.mol[j]
								uro.ifit_restart_needed(gfx) #Warning for ifit (restart automatically now)
								nmolformod -= 1
								break
					clean_mod_files(gfx,gfx.mod[i])
					Nma.remove_nmmod(gfx,gfx.mod[i])
					del gfx.mod[i]
					break
			gfx.renwin.Render()
		refresh_modlist(gfx,modlist)
		refresh_modlist(gfx,gfx.itf.modlist2)
		refresh_modlist(gfx,gfx.itf.modlist3)
		refresh_all_mollists(gfx)
		nfv.set('')
		modnameentry.delete(0,tk.END)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def refresh_modlist(gfx,mlst):
	try:
		mlst.delete(0, tk.END)
		for mod in gfx.mod:
			mlst.insert(tk.END,mod.un)
	except :
		pass
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def clean_mod_files(gfx,mod):
	chdir(gfx.workdir)
	if path.exists('d/xyz%d.d'%mod.id):
		system('rm -f d/xyz%d.d'%mod.id)
	if path.exists('d/map%d.d'%mod.id):
		system('rm -f d/map%d.d'%mod.id)
	if path.exists('o/tabl%d.s'%mod.id):
		system('rm -f o/tabl%d.s'%mod.id)
	if path.exists('f/tabl%d'%mod.id):
		system('rm -f f/tabl%d'%mod.id)
	if path.exists('f/itabl%d'%mod.id):
		system('rm -f f/itabl%d'%mod.id)
	if mod.rfn!=None:
		if path.exists(mod.rfn):
			system('rm -f %s'%mod.rfn)
	if mod.dfn!=None:
		if path.exists(mod.dfn):
			system('rm -f %s'%mod.dfn)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def set_mod_id(gfx):
	gfx.nmod += 1
	return gfx.nmod
#===============================================================================
#                              Defintion des molecule + wizard function
#===============================================================================
class Mol():
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def __init__(self):
		self.id=None
		self.fn=None
		self.un=None
		self.col=None
		self.mod=None
		self.acteur=None
		self.reader=None
		self.mapper=None
		self.rendtype=None
		self.bdsl=1.0 #bonds length
		self.hbdsl=1.0 #hydrogen bongs length
		self.lsm=[]#liste des symmate
		self.lnbsm=[]#list des numero de symetrie remplacer par l'intersection ou l'union de toutes
		self.slnbsm=[]#list specifique des numero de symetrie
		self.obs=None
		self.symobs=None
		self.nma=[]
		self.lclash = []
		self.lsphere=None
		self.iso=None
		self.nmcl=[]#Normal modes conformeres list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def which_mol(self,gfx,obj=None,event=None):
		wami=self.whoami(gfx)
		if gfx.whichmol!=wami:
			gfx.whichmol=wami
			gfx.itf.mollist.selection_clear(0,tk.END)
			gfx.itf.mollist.selection_set(wami-1)
			if gfx.ifit != None:
				gfx.ifit.change_mol()
		if gfx.ifit != None:
			gfx.ifit.fitin_io()
		if gfx.itf.mollist.curselection()==():
			gfx.itf.mollist.selection_clear(0,tk.END)
			gfx.itf.mollist.selection_set(wami-1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def whoami(self,gfx):
		for i in range(len(gfx.mol)):
			if gfx.mol[i].id == self.id:
				return i+1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def cords(self):
		mat=self.acteur.GetMatrix()
		(a,b,g) = R2Eul(ExtractRotMat(mat))
		(x,y,z) = [mat.GetElement(0,3),mat.GetElement(1,3),mat.GetElement(2,3)]
		return (a,b,g,x,y,z)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def move_lsphere(self,obj=None,event=None):
		self.lsphere.SetPosition(self.acteur.GetPosition())
		if self.render:
			self.gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Fonction hors class relative a la gestion des mols
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def chose_col(gfx,butcolor):
		new_color = tkColorChooser.askcolor (title="Molecule color")
		if new_color[1] != None:
				col = tk2vtk_color (new_color[0])
			    	butcolor.config(bg=vtk2tkhex_color(col))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def display_mol(gfx,mollist,butcolor):
	if mollist.curselection()==():
		MB.showwarning('Info','Select a molecule in the list')
		return
	for molit in mollist.curselection():
			for mol in gfx.mol:
				if mollist.get(molit)==mol.un:
					mol.col=tkhex2vtk_color(butcolor.cget('bg'))
					if mol.acteur.GetClassName()=='vtkAssembly':
						for i in range(mol.acteur.GetNumberOfPaths()):
							mol.acteur.GetParts().GetItemAsObject(i).GetProperty().SetColor(mol.col)
					else:
						mol.acteur.GetProperty().SetColor(mol.col)
					for sm in mol.lsm:
						if sm.GetClassName()=='vtkAssembly':
							for i in range(sm.GetNumberOfPaths()):
								sm.GetParts().GetItemAsObject(i).GetProperty().SetColor(Map.invcolor(mol.col))
					gfx.itf.mollist.itemconfigure(tk.END,background='black',foreground=vtk2tkhex_color(mol.col),selectbackground=vtk2tkhex_color(mol.col),selectforeground='white')
					mollist.itemconfigure(molit,background=butcolor.cget('bg'))
					gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def add_mol(gfx,mollist1,nameentry,modlist,extmollist,butcolor,rendvar,truecolor,varcolor,mollist2):
	molname=nameentry.get()
	if molname == '' :
		MB.showwarning('Info','Name entry empty, fill the blank')
		return
	mol=Mol()
	for molo in gfx.mol:
		if molo.un == molname:
			MB.showwarning('Info','Name already defined, chose another one')
			return
	mol.un=molname
	for selmod in modlist.curselection():
			for mod in gfx.mod:
				if modlist.get(selmod)==mod.un:
					mol.mod=mod
	try :
		mol.mod.fn
	except AttributeError:
		MB.showwarning('Info','Select a model in the list. If the list is empty go to File|Models ')
		return
	mol.id=set_mol_id(gfx)
	mol.mod.nmol +=1
	gfx.whichmol=mol.id
	mol.col=tkhex2vtk_color(butcolor.cget('bg'))
	mollist1.insert(tk.END, mol.un)
	mollist1.itemconfigure(tk.END,background=butcolor.cget('bg'))
	gfx.mol+=[mol]
	rendic={'Line':0, 'Ball':1,'Stick':2,'Ball & Line':4,'Ball & Stick':3,'Surface':5,'Wireframe':6,'Points':7}
	mol.rendtype=rendic[rendvar.get()]
	gfx.itf.root.configure(cursor='watch')
	gfx.itf.status.set('Molecule loading ... please wait')
	if mol.mod.type=='mol':
		if mol.mod.dspmodtype == 'CA':
			mol.bdsl=4.0
			mol.hdbsl=0.0
			load_mols(mol,rendic[rendvar.get()],gfx,truecolor.get())
		elif mol.mod.dspmodtype == 'AA' or mol.mod.dspmodtype == 'B':
			mol.bdsl=1.0
			mol.hdbsl=1.0
			load_mols(mol,rendic[rendvar.get()],gfx,truecolor.get())
		else :
			print 'strange Error in mol.mod.dspmodtype'
	elif mol.mod.type=='map':
		load_maps(mol,rendic[rendvar.get()],gfx,truecolor.get())
	else :
		print 'strange Error in mol.mod.type'
	mol.obs= mol.acteur.AddObserver("ModifiedEvent",lambda obj,event,mol=mol,gfx=gfx : mol.which_mol(gfx,obj,event))
	moddic={}
	for mol in gfx.mol:
		if mol.mod.un in moddic.keys():
			moddic[mol.mod.un]+=1
		else:
			moddic[mol.mod.un]=1
	if modlist.get(tk.ACTIVE) not in moddic.keys():
		nameentry.delete(0,tk.END)
		nameentry.insert(tk.END,'%s_1'%modlist.get(tk.ACTIVE))
	else :
		nameentry.delete(0,tk.END)
		nameentry.insert(tk.END,'%s_%d'%(modlist.get(tk.ACTIVE),(moddic[modlist.get(tk.ACTIVE)]+1)))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def set_om_color(item,col,tc):
			varcolor.set(item)
			butcolor.configure(background=col)
			truecolor.set(tc)
	allcol =[('red','#ff393c',0),('green','#46ff76',0),('purple','#ff5cec',0),('yellow','#f9f341',0),('indigo','#5cf3e6',0),('blue','#0098d8',0)]
	numcol = gfx.nmol%6
	pcol = allcol[numcol]
	set_om_color(pcol[0],pcol[1],pcol[2])
	refresh_all_mollists(gfx)
	gfx.itf.status.clear()
	gfx.itf.root.configure(cursor='arrow')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def refresh_mollist(gfx,mlst,ext):
		mlst.delete(0, tk.END)
		for mol in gfx.mol:
			mlst.insert(tk.END,mol.un)
			if ext == 1:
				mlst.itemconfigure(tk.END,background='black',foreground=vtk2tkhex_color(mol.col),selectbackground=vtk2tkhex_color(mol.col),selectforeground='white')
			else :
				mlst.itemconfigure(tk.END,background=vtk2tkhex_color(mol.col))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def refresh_all_mollists(gfx):
	loml=[gfx.itf.mollist2,gfx.itf.mollist3,gfx.itf.mollist4]
	refresh_mollist(gfx,gfx.itf.mollist,1)
	Nma.build_lists(gfx,gfx.itf.nmmodlist,gfx.itf.nmmollist)
	for mlst in loml:
		try:
			refresh_mollist(gfx,mlst,0)
		except:
			pass
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def remove_mol(gfx,mollist1,extmollist,mollist2):
	if gfx.mol != []:
		if mollist1.curselection() == ():
			MB.showwarning('Info','Select a molecule in the list')
		for selmol in mollist1.curselection():
			for i in range(len(gfx.mol)):
				if mollist1.get(selmol)==gfx.mol[i].un:
					for sm in gfx.mol[i].lsm:
						gfx.renderer.RemoveActor(sm)
						del sm
					gfx.renderer.RemoveActor(gfx.mol[i].acteur)
					if gfx.mol[i].lsphere!=None:
						gfx.renderer.RemoveActor(gfx.mol[i].lsphere)
					del gfx.mol[i]
					uro.ifit_restart_needed(gfx) #Warning for ifit
					gfx.renwin.Render()
					break
		refresh_all_mollists(gfx)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def up_or_down_mol(gfx,mollist1,extmollist,mollist2,sens):#sens = 1 for up | -1 for down
	if gfx.mol != []:
		if mollist1.curselection() == ():
			MB.showwarning('Info','Select a molecule in the list')
		selection =  list(mollist1.curselection())
		if sens==-1:
			selection.reverse()
		for selmol in selection:
			for i in range(len(gfx.mol)):
				if mollist1.get(selmol)==gfx.mol[i].un:
					tmp = gfx.mol[(i-1*sens)%len(gfx.mol)]
					gfx.mol[(i-1*sens)%len(gfx.mol)]=gfx.mol[i]
					gfx.mol[i]=tmp
					break
		if selection != []:
			refresh_all_mollists(gfx)
			uro.ifit_restart_needed(gfx) #Warning for ifit
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def set_mol_id(gfx):
	gfx.nmol += 1
	return gfx.nmol
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def printmolpv(gfx):
	for mol in gfx.mol:
		mat=mol.acteur.GetMatrix()
		(a,b,g) = sym.R2Eul(sym.ExtractRotMat(mat))
		(x,y,z) = [mat.GetElement(0,3),mat.GetElement(1,3),mat.GetElement(2,3)]
		print' mol : %d : %.3f %.3f %.3f %.3f %.3f %.3f\n'%(mol.id,a,b,g,x,y,z)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def hide_mol(gfx,mollist):
	if mollist.curselection()==():
		MB.showwarning('Info','Select a molecule in the list')
		return
	for molit in mollist.curselection():
			for mol in gfx.mol:
				if mollist.get(molit)==mol.un:
					mol.acteur.VisibilityOff()
	gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def show_mol(gfx,mollist):
	if mollist.curselection()==():
		MB.showwarning('Info','Select a molecule in the list')
		return
	for molit in mollist.curselection():
			for mol in gfx.mol:
				if mollist.get(molit)==mol.un:
					mol.acteur.VisibilityOn()
	gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def hide_constellation(gfx,mollist):
	if mollist.curselection()==():
		MB.showwarning('Info','Select a molecule in the list')
		return
	for molit in mollist.curselection():
			for mol in gfx.mol:
				if mollist.get(molit)==mol.un:
					if mol.lsm!=[]:
						for sm in mol.lsm:
							sm.VisibilityOff()
						gfx.renwin.Render()
					else :
						MB.showwarning('Info','Missing constellation for one molecule')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def show_constellation(gfx,mollist):
	if mollist.curselection()==():
		MB.showwarning('Info','Select a molecule in the list')
		return
	for molit in mollist.curselection():
			for mol in gfx.mol:
				if mollist.get(molit)==mol.un:
					if mol.lsm!=[]:
						for sm in mol.lsm:
							sm.VisibilityOn()
						gfx.renwin.Render()
					else :
						MB.showwarning('Info','Missing constellation for one molecule')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def invert_constellation(gfx):
	try:
		if gfx.mol[0].lsm[0].GetVisibility():
			for mol in gfx.mol:
				for sm in mol.lsm:
					sm.VisibilityOff()
		else:
			for mol in gfx.mol:
				for sm in mol.lsm:
					sm.VisibilityOn()
		gfx.renwin.Render()
	except :
		pass
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def locate_independent_mol(gfx):
	mode_localisation=0
	if gfx.mol==[]:
		return
	for mol in gfx.mol:
		if mol.lsphere!=None:
			mode_localisation=1
	if mode_localisation==0:
		source = vtk.vtkSphereSource()
		source.SetThetaResolution(8)
		source.SetPhiResolution(8)
		source.SetRadius(mol.mod.sphdist)
		mapper = vtk.vtkPolyDataMapper()
		mapper.SetInput(source.GetOutput())
		for mol in gfx.mol:
			sphere = vtk.vtkActor()
			sphere.SetMapper(mapper)
			sphere.GetProperty().SetColor(1,1,1)
			sphere.GetProperty().SetOpacity(0.5)
			sphere.PickableOff()
			sphere.DragableOff()
			sphere.SetUserMatrix(mol.acteur.GetMatrix())
			mol.lsphere = sphere
			gfx.renderer.AddActor(mol.lsphere)
	else :
		for mol in gfx.mol:
			gfx.renderer.RemoveActor(mol.lsphere)
			del mol.lsphere
			mol.lsphere=None
	gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def build_bookkeeping(gfx):
	bk = open('%s/mols.bk'%gfx.tmpdir,'w')
	bk.write('%3s  %7s%20s%12s%20s     %-50s%-s\n'%('#','Mol-id','Mol-name','Mod-id','Mod-name','Model-filename in ./%s/'%Map.extract_file_from_path(gfx.tmpdir),'Source-model-filename'))
	i=0
	for mol in gfx.mol:
		i+=1
		bk.write('%3s  %7s%20s%12s%20s     %-50s%-s\n'%(i,str(mol.id),mol.un, str(mol.mod.id),mol.mod.un,Map.extract_file_from_path(mol.mod.rfn),mol.mod.ofn))
	bk.close()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def isca(pdbfile):
	pdb = open(pdbfile,'r')
	nbcon=0.
	nbca=0.
        for l in pdb:
		if 'ATOM' in l:
		 	if ' CA  ' in l :
		        	nbca+=1
			elif  ' C  ' in l or ' N  ' in l or ' O  ' in l  :
				nbcon+=1
	pdb.close()
	if nbcon==0 :
		return True
	else :
		return nbca/nbcon > 0.75
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              LOAD ALL ATOM and CA SUB UNIT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_mols(mol,rendmod,gfx,atomcol):
	mol.reader = vtk.vtkPDBReader()
	mol.reader.SetFileName(mol.mod.dfn)
	mol.reader.Update() #by calling Update() we read the file
	mol.reader.SetHBScale(mol.hbdsl)
	mol.reader.SetBScale(mol.bdsl)
	Localmol=[]
	beg=0
	if mol.acteur!=None:
		pos = mol.acteur.GetPosition()
		ori = mol.acteur.GetOrientation()
		gfx.renderer.RemoveActor(mol.acteur)
		beg=1
	if rendmod in [0,4] :
		mol.mapper = vtk.vtkPolyDataMapper()
		mol.mapper.SetInputConnection(mol.reader.GetOutputPort())
		if atomcol!=1:
			mol.mapper.UseLookupTableScalarRangeOff()
			mol.mapper.SetScalarVisibility(0)
		else:
			mol.mapper.SetScalarModeToDefault()
		lacteur = vtk.vtkActor()
		lacteur.GetProperty().SetColor(mol.col)
		if mol.mod.dspmodtype=='AA':
			lacteur.GetProperty().SetLineWidth(3)
		else :
			lacteur.GetProperty().SetLineWidth(4)
		lacteur.SetMapper(mol.mapper)
		Localmol+=[lacteur]
	if rendmod in [1,3,4]:
		sphere = vtk.vtkSphereSource()
		sphere.SetCenter(0, 0, 0)
		sphere.SetRadius(1)
		sphere.SetThetaResolution(8)
		sphere.SetStartTheta(0)
		sphere.SetEndTheta(360)
		sphere.SetPhiResolution(8)
		sphere.SetStartPhi(0)
		sphere.SetEndPhi(180)
		glyph = vtk.vtkGlyph3D()
		glyph.SetInputConnection(mol.reader.GetOutputPort())
		glyph.SetColorMode(1)
		glyph.ScalingOn()
		glyph.SetScaleMode(2)
		glyph.SetScaleFactor(0.25)
		glyph.SetSource(sphere.GetOutput())
		mol.mapper = vtk.vtkPolyDataMapper()
		mol.mapper.SetInputConnection(glyph.GetOutputPort())
		if atomcol!=1:
			mol.mapper.UseLookupTableScalarRangeOff()
			mol.mapper.SetScalarVisibility(0)
		else:
			mol.mapper.SetScalarModeToDefault()
		bacteur = vtk.vtkActor()
		bacteur.SetMapper(mol.mapper)
		bacteur.GetProperty().SetRepresentationToSurface()
		bacteur.GetProperty().SetInterpolationToGouraud()
		bacteur.GetProperty().SetAmbient(0.15)
		bacteur.GetProperty().SetDiffuse(0.85)
		bacteur.GetProperty().SetSpecular(0.1)
		bacteur.GetProperty().SetSpecularPower(100)
		bacteur.GetProperty().SetSpecularColor(1, 1, 1)
		bacteur.GetProperty().SetColor(mol.col)
		Localmol+=[bacteur]
	if rendmod in [2,3] :
		tubes = vtk.vtkTubeFilter()
		tubes.SetInputConnection(mol.reader.GetOutputPort())
		tubes.SetNumberOfSides(8)
		tubes.SetCapping(0)
		tubes.SetRadius(1)
		tubes.SetVaryRadius(0)
		tubes.SetRadiusFactor(10)
		mol.mapper = vtk.vtkPolyDataMapper()
		mol.mapper.SetInputConnection(tubes.GetOutputPort())
		if atomcol!=1:
			mol.mapper.UseLookupTableScalarRangeOff()
			mol.mapper.SetScalarVisibility(0)
		else:
			mol.mapper.SetScalarModeToDefault()
		tacteur = vtk.vtkActor()
		tacteur.SetMapper(mol.mapper)
		tacteur.GetProperty().SetRepresentationToSurface()
		tacteur.GetProperty().SetInterpolationToGouraud()
		tacteur.GetProperty().SetAmbient(0.15)
		tacteur.GetProperty().SetDiffuse(0.85)
		tacteur.GetProperty().SetSpecular(0.1)
		tacteur.GetProperty().SetSpecularPower(100)
		tacteur.GetProperty().SetSpecularColor(1, 1, 1)
		tacteur.GetProperty().SetColor(mol.col)
		Localmol+=[tacteur]
	assembly = vtk.vtkAssembly()
	for act in Localmol:
		assembly.AddPart(act)
		del act
	mol.acteur=assembly
	if beg:
		mol.acteur.SetPosition(pos)
		mol.acteur.SetOrientation(ori)
	else :
		remember_rottra(mol)
	gfx.renderer.AddActor(mol.acteur)
	gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def remember_rottra(mol):
	a,b,g = mol.mod.rottra[0],mol.mod.rottra[1],mol.mod.rottra[2]
	x,y,z = mol.mod.rottra[3],mol.mod.rottra[4],mol.mod.rottra[5]
	RotaEuler(mol.acteur,a,b,g)
	mol.acteur.SetPosition(x,y,z)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def RotaEuler(acteur,alpha,beta,gamma):
		acteur.RotateWXYZ(gamma,0,0,1)
		acteur.RotateWXYZ(beta,0,1,0)
		acteur.RotateWXYZ(alpha,0,0,1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              LOAD Maps molecules
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_maps(mol,rendmod,gfx,atomcol):
	mol.reader = vtk.vtkStructuredPointsReader()
	mol.reader.SetFileName(mol.mod.dfn)
	mol.reader.Update() #by calling Update() we read the file
	mol.iso = vtk.vtkMarchingContourFilter()
	mol.iso.UseScalarTreeOn()
	mol.iso.ComputeNormalsOn()
	mol.iso.SetInputConnection(mol.reader.GetOutputPort())
	mol.iso.SetValue(0,mol.mod.isov*mol.mod.sigavg[0]+mol.mod.sigavg[1])
	clean = vtk.vtkCleanPolyData()
  	clean.SetInputConnection(mol.iso.GetOutputPort())
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
	mol.mapper = vtk.vtkOpenGLPolyDataMapper()
	mol.mapper.SetInputConnection(smooth.GetOutputPort()) ### <- connection here
	mol.mapper.ScalarVisibilityOff()
	mol.mapper.Update()
	mol.acteur= vtk.vtkOpenGLActor()
	mol.acteur.SetMapper(mol.mapper)
	mol.acteur.GetProperty().SetColor(mol.col)
	if rendmod==5:
		mol.acteur.GetProperty().SetRepresentationToSurface()
	elif rendmod==6:
		mol.acteur.GetProperty().SetRepresentationToWireframe()
	elif rendmod==7:
		mol.acteur.GetProperty().SetRepresentationToPoints()
	else :
		mol.acteur.GetProperty().SetRepresentationToSurface()
	mol.acteur.GetProperty().SetInterpolationToGouraud()
	mol.acteur.GetProperty().SetSpecular(.4)
	mol.acteur.GetProperty().SetSpecularPower(10)
	gfx.renderer.AddActor(mol.acteur)
	gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              conversion de couleur
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def tk2vtk_color(tk_col):
	""" conversion of color codes from tk to vtk"""
	ONE_255 = 1.0/255.0
	return (tk_col[0]*ONE_255, tk_col[1]*ONE_255, tk_col[2]*ONE_255)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def vtk2tk_color (vtk_col):
	""" conversion of color codes from vtk to tk"""
	return (vtk_col[0]*255.0, vtk_col[1]*255.0, vtk_col[2]*255.0)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def vtk2tkhex_color(vtk_col):
	(r,v,b)=(vtk_col[0]*255.0, vtk_col[1]*255.0, vtk_col[2]*255.0)
	tk_rgb = "#%02x%02x%02x" % (r, v, b)
  	return tk_rgb
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def tkhex2vtk_color(tk_hex_col):
	split = (tk_hex_col[1:3], tk_hex_col[3:5], tk_hex_col[5:7])
	(r,g,b)= [int("0x"+x, 16) for x in split]
	return (r/255.,g/255.,b/255.)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Save molecule as
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def save_mol_as(root,mollist,gfx):
	if mollist.curselection()==():
		MB.showwarning('Info','Select a molecule in Mol List for saving it')
		return
	for molit in mollist.curselection():
		for mol in gfx.mol:
			if mollist.get(molit)==mol.un:
				fn = FD.asksaveasfilename(parent=root,filetypes=[('PDB files','*.pdb')] ,title="Save molecule %s coordinates as..."%mol.un)
				if len(fn ) > 0:
					if not fn.endswith('.pdb'):
						fn = fn +'.pdb'
		    			mat = mol.acteur.GetMatrix() #recuperation de la mat de pv 4x4
		    			R = array(sym.ExtractRotMat(mat)) #extraction de la matrice de rotation depuit la vtk 4x4
					t = array([mat.GetElement(0,3),mat.GetElement(1,3),mat.GetElement(2,3)]) #extract tx,ty,tz
					pdbin=open(mol.mod.fn,"r")
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
					#ferme les fichiers
					pdbin.close()
					pdbout.close()
					return
#===============================================================================
#                              mini renderer windows for models
#===============================================================================
class miniwin():
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def __init__(self,frame):
		self.renderer = vtk.vtkOpenGLRenderer()
		self.renwin = vtk.vtkXOpenGLRenderWindow()
		self.renwin.AddRenderer(self.renderer)
		self.renwin.SetStereoTypeToAnaglyph()
		#self.renwin.SetAnaglyphColorMask(4,2)
		self.iren = vtkTkRenderWindowInteractor(frame,rw=self.renwin, width=300, height=300 )
		self.istyle = vtk.vtkInteractorStyleSwitch()
		self.iren.SetInteractorStyle(self.istyle)
		self.istyle.SetCurrentStyleToTrackballCamera()
		self.iren.Initialize()
		self.iren.pack(side='bottom', fill='none', expand=0)
		self.iren.Start()
		self.camera=vtk.vtkCamera()
		self.camera.SetFocalPoint(0, 0, 0)
		self.camera.SetPosition(0, 0, 250)
		self.camera.SetViewUp(0, 0, 0)
		self.camera.SetEyeAngle(5.0)
		self.renderer.SetActiveCamera(self.camera)
		self.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def display_rottra(self,mod):
		tm = vtk.vtkTextMapper()
		tp = tm.GetTextProperty()
		#tp.SetFontFamilyToArial()
		tp.SetFontSize(10)
		#tp.BoldOff()
		tp.ShadowOff()
		tp.SetColor(1, 1, 1)
		tp.SetOpacity(1)
		label = vtk.vtkActor2D()
		label.VisibilityOn()
		label.SetMapper(tm)
		v=[0,15]
		label.SetPosition(v)
		rottra = mod.rottra
		tm.SetInput("R:[%.3f,%.3f,%.3f]  T:[%.3f,%.3f,%.3f]"%(rottra[0],rottra[1],rottra[2],rottra[3],rottra[4],rottra[5]))
		self.renderer.AddActor2D(label)
		self.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def display_mod(self,mod):
		self.renderer.RemoveAllViewProps()
		if mod.type == 'mol':
			if mod.dspmodtype == 'CA':
				reader = vtk.vtkPDBReader()
				reader.SetFileName(mod.dfn)
				reader.SetHBScale(0)
				reader.SetBScale(4.0)
				mapper = vtk.vtkPolyDataMapper()
				mapper.SetInputConnection(reader.GetOutputPort())
				mapper.SetScalarModeToDefault()
				acteur = vtk.vtkActor()
				acteur.SetMapper(mapper)
			else: #all atoms or backbone case
				reader = vtk.vtkPDBReader()
				reader.SetFileName(mod.dfn)
				mapper = vtk.vtkPolyDataMapper()
				mapper.SetInputConnection(reader.GetOutputPort())
				mapper.SetScalarModeToDefault()
				acteur = vtk.vtkActor()
				acteur.SetMapper(mapper)
			self.camera.SetFocalPoint(0, 0, 0)
			self.camera.SetPosition(0, 0, 250)
		elif mod.type=='map':
			reader = vtk.vtkStructuredPointsReader()
			reader.SetFileName(mod.dfn)
			reader.Update() #by calling Update() we read the file
			iso = vtk.vtkMarchingContourFilter()
			iso.UseScalarTreeOn()
			iso.ComputeNormalsOn()
			iso.SetInputConnection(reader.GetOutputPort())
			iso.SetValue(0,mod.isov*mod.sigavg[0]+mod.sigavg[1])
			mapper = vtk.vtkOpenGLPolyDataMapper()
			mapper.SetInputConnection(iso.GetOutputPort()) ### <- connection here
			mapper.ScalarVisibilityOff()
			mapper.Update()
			acteur= vtk.vtkOpenGLActor()
			acteur.SetMapper(mapper)
			acteur.GetProperty().SetColor(0,0.35,1)
			if mod.rep=='Wireframe':
				acteur.GetProperty().SetRepresentationToWireframe()
			elif mod.rep=='Surface':
				acteur.GetProperty().SetRepresentationToSurface()
			else :
				acteur.GetProperty().SetRepresentationToWireframe()
			acteur.GetProperty().SetInterpolationToGouraud()
			acteur.GetProperty().SetSpecular(.4)
			acteur.GetProperty().SetSpecularPower(10)
			(xmin,xmax,ymin,ymax,zmin,zmax)=acteur.GetBounds()
			maxi=max(xmax-xmin,ymax-ymin,zmax-zmin)
			mini=min(xmax-xmin,ymax-ymin,zmax-zmin)
			fp=((xmax+xmin)/2.,(ymax+ymin)/2.,(zmax+zmin)/2.)
			self.camera.SetFocalPoint(fp[0],fp[1],fp[2])
			self.camera.SetPosition(fp[0],fp[1],fp[2]+(zmax-zmin)*4.)
		else :
			print 'strange Error in mod.type'
		(xl, xu, yl, yu, zl, zu)=acteur.GetBounds()
		mod.coldist = min(xu-xl,yu-yl,zu-zl)/2.
		mod.sphdist = max(xu-xl,yu-yl,zu-zl)/2.
		self.renderer.AddActor(acteur)
		self.renderer.ResetCameraClippingRange()
		self.renwin.Render()
