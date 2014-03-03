# -*- coding: utf-8 -*-
# Symmetry Module of VEDA
#
# Copyright:
# 2009-01-01 Gael Goret  gael.goret@ibs.fr
#
# Last modified:
# 2010-09-06 Gael Goret  gael.goret@ibs.fr
try :
	import Tkinter as tk
except :
	import tkinter as tk
import vtk
import gfx,Map,mod,itf,uro
from numpy import *
from numpy.linalg import det
from math import *
from os import chdir, environ, system, mkdir
import tkMessageBox as MB
import tkFileDialog as FD
#===============================================================================
class symmate(vtk.vtkAssembly):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def __init__(self):
		self.ut = vtk.vtkTransform()
		self.ut.PostMultiply()
		self.nbsym=None
#===============================================================================
class Sym():
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def __init__(self,mollist,gfx,mdbe,indi,coli):
		self.mollist=mollist
		self.symlistfile=gfx.workdir + '/SYM'
		self.gfx=gfx
		self.mdbe=mdbe #Minimal distance to box edges
		self.symlist=[]#intersec symlist
		self.unionsymlist=[]
		self.render = 1
		self.indi = int(indi)
		self.coli = int(coli)
		if gfx.map!=[]:
			self.bounds=gfx.map[0].box.GetBounds()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def build_mol_symlist(self):
		for mol in self.gfx.mol:
			fst=1
			mol.lnbsm=[]
			mol.slnbsm=[]
			mol.lclash=[]
			rips=[] #real instance position vector
			(xmin, xmax, ymin, ymax, zmin, zmax)= self.bounds
			sym=open(self.symlistfile,'r')
			for l in sym:
				ms = l.split()
				nbl = int(ms[6][1:])
				ang = [float(ms[0]),float(ms[1]),float(ms[2])]
				tra = array([float(ms[3]),float(ms[4]),float(ms[5])])
				sm=symmate()
				sm.SetPosition(mol.acteur.GetPosition())
				sm.SetOrientation(mol.acteur.GetOrientation())
				self.RotaEuler(sm.ut,ang[0],ang[1],ang[2])
				sm.ut.Translate(tra[0],tra[1],tra[2])
				sm.SetUserTransform(sm.ut)
				pip = [sm.GetMatrix().GetElement(0,3),sm.GetMatrix().GetElement(1,3),sm.GetMatrix().GetElement(2,3)]
				if (xmin + self.mdbe < pip[0]) and (pip[0] < xmax - self.mdbe) and (ymin + self.mdbe < pip[1]) and (pip[1] < ymax - self.mdbe) and (zmin + self.mdbe < pip[2]) and (pip[2] < zmax - self.mdbe):
					mol.lnbsm+=[nbl]
					mol.slnbsm+=[nbl]
					collision=0
					for rip in rips:
						if sqrt( (pip[0]-rip[0])**2 + (pip[1]-rip[1])**2 + (pip[2]-rip[2])**2) < mol.mod.coldist:
							mol.lclash+=[nbl]
							collision=1
							break
					if not collision:
						rips+=[pip]
			sym.close()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def remove_mol_clash(self):
		for mol in self.gfx.mol:
			for c in mol.lclash:
				try: # si non self.indi c's peut ne pas etre present
					mol.lnbsm.pop(mol.lnbsm.index(c))
					mol.slnbsm.pop(mol.slnbsm.index(c))
				except: pass
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def build_intersect_symlist(self):
		symlist=[]
		for mol in self.gfx.mol:
			if symlist == []:
				symlist = mol.slnbsm[:]
				continue
			intersec=[]
			for i in mol.slnbsm:
				if i in symlist:
					intersec += [i]
			symlist=intersec
		self.symlist=symlist
		for mol in self.gfx.mol:
			if self.indi :
				mol.lnbsm=mol.slnbsm[:]
			else :
				mol.lnbsm=self.symlist[:]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def build_union_symlist(self):
		union=[]
		for mol in self.gfx.mol:
			if union == []:
				union=mol.slnbsm[:]
				continue
			for i in mol.slnbsm:
				if i not in union:
					union += [i]
		union.sort()
		self.unionsymlist=union
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def build_symlist_file(self):
		syml = open('symlist','w')
		if self.indi:
			smlst = ' '+ str(self.unionsymlist).replace(',','')[1:-1]+'\n'
		else :
			smlst = ' '+ str(self.symlist).replace(',','')[1:-1]+'\n'
		syml.write(smlst)
		syml.close()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              building
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def build_constelation(self,mol):
		if mol.acteur == None :
			MB.showwarning('Info','Select a molecule in the list')
			return
		if mol.symobs !=None:
			mol.acteur.RemoveObserver(mol.symobs)
		if mol.lsm!=[]:
			for sm in mol.lsm:
				self.gfx.renderer.RemoveActor(sm)
			mol.lsm=[]
		(xmin, xmax, ymin, ymax, zmin, zmax)= self.bounds
		sym=open(self.symlistfile,'r')
		for l in sym:
			ms = l.split()
			nbl = int(ms[6][1:])
			if nbl not in mol.lnbsm:
				continue
			ang = [float(ms[0]),float(ms[1]),float(ms[2])]
			tra = array([float(ms[3]),float(ms[4]),float(ms[5])])
			sm=symmate() #on cree un symmate vide
			sm.SetPosition(mol.acteur.GetPosition()) #on assigne la partie translationelle des pv
			sm.SetOrientation(mol.acteur.GetOrientation()) #on assigne la partie rotationelle des pv
			self.RotaEuler(sm.ut,ang[0],ang[1],ang[2]) #on defini la partie rotationelle de la transformation
			sm.ut.Translate(tra[0],tra[1],tra[2]) #on defini la partie translationelle de la transformation
			sm.SetUserTransform(sm.ut) #on assigne la transformation a notre symmate
			pip = [sm.GetMatrix().GetElement(0,3),sm.GetMatrix().GetElement(1,3),sm.GetMatrix().GetElement(2,3)]#on recupere la partie translationelle de la combinaison de pv et de la transformation (ut)
			if (xmin + self.mdbe < pip[0]) and (pip[0] < xmax - self.mdbe) and (ymin + self.mdbe < pip[1]) and (pip[1] < ymax - self.mdbe) and (zmin + self.mdbe < pip[2]) and (pip[2] < zmax - self.mdbe):# on test si pip est dans la boite
				sm.nbsym=nbl
				if mol.acteur.GetClassName()=='vtkAssembly':# dans le cas ou la molecule independante est un assembly
					for i in range(mol.acteur.GetNumberOfPaths()):
						tmp=vtk.vtkActor()
						tmp.SetMapper(mol.acteur.GetParts().GetItemAsObject(i).GetMapper())
						p=vtk.vtkProperty()
						#p.SetColor(mol.acteur.GetParts().GetItemAsObject(i).GetProperty().GetColor())
						p.SetColor(Map.invcolor(mol.acteur.GetParts().GetItemAsObject(i).GetProperty().GetColor()))
						tmp.SetProperty(p)
						if mol.mod.type=='mol':
							tmp.GetProperty().SetLineWidth(4)
						tmp.DragableOff()
						tmp.PickableOff()
						sm.AddPart(tmp)
				else:#cas simple ou la mol ind est composer d un seul objet
					tmp=vtk.vtkActor()
					tmp.SetMapper(mol.acteur.GetMapper())
					p=vtk.vtkProperty()
					#p.SetColor(mol.acteur.GetParts().GetItemAsObject(i).GetProperty().GetColor())
					p.SetColor(Map.invcolor(mol.acteur.GetProperty().GetColor()))
					tmp.SetProperty(p)
					if mol.mod.type=='mol':
							tmp.GetProperty().SetLineWidth(4)
					tmp.DragableOff()
					tmp.PickableOff()
					sm.AddPart(tmp)
				mol.lsm+=[sm]# on ajoute le symmate a la liste des symmate
		sym.close()
		self.move_sym(mol)
#                              euler rotation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def RotaEuler(self,acteur,alpha,beta,gamma):
		acteur.RotateWXYZ(gamma,0,0,1)
		acteur.RotateWXYZ(beta,0,1,0)
		acteur.RotateWXYZ(alpha,0,0,1)
#                              render
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def render_constelation(self,mol):
		if mol.lsm!=[]:
			for sm in mol.lsm:
				self.gfx.renderer.AddActor(sm)
			self.gfx.renderer.RemoveActor(mol.acteur)
			self.gfx.renderer.AddActor(mol.acteur)
			self.gfx.renwin.Render()
#                              event
	#definissons les fonctions qui feront bouger les shallowcopy en mm temps que la molecule independante :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def add_event(self,mol):
		mol.symobs= mol.acteur.AddObserver('ModifiedEvent',lambda obj,event,mol=mol :self.move_sym(mol,obj,event))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def move_sym(self,mol,obj=None,event=None):
		for sm in mol.lsm :
			sm.SetPosition(mol.acteur.GetPosition())
			sm.SetOrientation(mol.acteur.GetOrientation())
		if self.render:
			self.gfx.renwin.Render()
#                              refresh fonction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def refresh_all_constelation(self,gfx):
		gfx.itf.root.configure(cursor='watch')
		gfx.itf.status.set('Constellation loading ... please wait')
		if gfx.map!=[]:
			self.bounds=gfx.map[0].box.GetBounds() #update des bounds
		self.build_mol_symlist() #general cree tout les sym list de toutes les mol
		self.build_intersect_symlist() #calcule l'intersection de tout les sym -> self.symlist
		self.build_union_symlist() #calcule l'union de tout les sym -> self.unionsymlist
		if self.coli:
			self.remove_mol_clash() #enleve de la list des molecule symetrique celle qui sont en clash
		self.build_symlist_file() #construit le fichier symlist
		if gfx.ifit != None:
			uro.reset_ifit(gfx)
		for mol in gfx.mol:
			self.build_constelation(mol) #on cree la constelation (la liste de symmate)
			self.render_constelation(mol) #on afficher tout ce petit monde
			self.add_event(mol) #on ajoute l event links
		gfx.itf.status.clear()
		gfx.itf.root.configure(cursor='arrow')
#                              MAIN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def main(self):
		if self.mollist.curselection()==():
			MB.showwarning('Info','Select a molecule in the list')
			return
		if self.gfx.map==[]:
			MB.showwarning('Info','Load a map')
			return
		if self.gfx.ps==None:
			MB.showwarning('Info','Define symmetry')
			return
		self.gfx.itf.root.configure(cursor='watch')
		self.gfx.itf.status.set('Constellation loading ... please wait')
		self.build_mol_symlist() #general cree tout les sym list de toutes les mol
		self.build_intersect_symlist() #calcule l'intersection de tout les sym -> self.symlist
		self.build_union_symlist() #calcule l'union de tout les sym -> self.unionsymlist
		if self.coli:
			self.remove_mol_clash() #enleve de la list des molecule symetrique celle qui sont en clash
		self.build_symlist_file() #construit le fichier symlist
		if self.gfx.ifit != None:
			uro.reset_ifit(self.gfx)
		for molit in self.mollist.curselection():
			for mol in self.gfx.mol:
				if self.mollist.get(molit)==mol.un:
					self.build_constelation(mol) #on cree la constelation (la liste de symmate)
					self.render_constelation(mol) #on afficher tout ce petit monde
					self.add_event(mol) #on ajoute l event links
		self.gfx.itf.status.clear()
		self.gfx.itf.root.configure(cursor='arrow')
#                              SAVE CONSTELLATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def save_constellation_as(self,mollist):
		if mollist.curselection()==():
			MB.showwarning('Info','Select a molecule in Mol List for saving it')
			return
		for molit in mollist.curselection():
			for mol in self.gfx.mol:
				if mollist.get(molit)==mol.un:
					fn = FD.asksaveasfilename(parent=self.gfx.itf.root,filetypes=[('PDB files','*.pdb')] ,title="Save constellation of molecule %s as..."%mol.un)
					if len(fn ) > 0:
						self.gfx.itf.root.configure(cursor='watch')
						self.gfx.itf.status.set('Save Constellation in progress ... please wait')
						if not fn.endswith('.pdb'):
							fn = fn +'.pdb'
						foldname = self.extract_file_from_path(fn)[:-4]
						foldpath = self.extract_path(fn)
						fold = foldpath + foldname
						try : mkdir(fold)
						except OSError :
							if not MB.askyesno('Info','Overwrite directory %s ?'%foldname):
								return
						pdbout=open(fn,"w")
						chdir(fold)
						mat = mol.acteur.GetMatrix()
			    			R = array(ExtractRotMat(mat))
						t = array([mat.GetElement(0,3),mat.GetElement(1,3),mat.GetElement(2,3)])
						pdbin=open(mol.mod.fn,"r")
						pdbconst=open('%s/1_%s.pdb'%(fold,foldname),"w")
						for ligne in pdbin:
							if (ligne[0:6] == "ATOM  " or ligne[0:6] == "HETATM"):
								x = float(ligne[30:38])
								y = float(ligne[38:46])
								z = float(ligne[46:54])
								r = array([x,y,z])
								newr = dot(R,r) + t
								ligne=ligne[:30]+"%8.3f%8.3f%8.3f"%tuple(newr)+ligne[54:]
							pdbconst.write(ligne)
						pdbconst.close()
						for sm in mol.lsm:
				    			mat = sm.GetMatrix()
				    			R = array(ExtractRotMat(mat))
							t = array([mat.GetElement(0,3),mat.GetElement(1,3),mat.GetElement(2,3)])
							pdbconst=open('%s/%d_%s.pdb'%(fold,sm.nbsym,foldname),"w")
							pdbin=open(mol.mod.fn,"r")
							for ligne in pdbin:
								if (ligne[0:6] == "ATOM  " or ligne[0:6] == "HETATM"):
									x = float(ligne[30:38])
									y = float(ligne[38:46])
									z = float(ligne[46:54])
									r = array([x,y,z])
									newr = dot(R,r) + t
									ligne=ligne[:30]+"%8.3f%8.3f%8.3f"%tuple(newr)+ligne[54:]
								pdbconst.write(ligne)
								pdbout.write(ligne)
							pdbconst.close()
							pdbin.close()
						#ferme les fichiers
						pdbout.close()
						chdir(self.gfx.workdir)
						self.gfx.itf.status.clear()
						self.gfx.itf.root.configure(cursor='arrow')
						return
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def extract_file_from_path(self,fn):
		for i in range(len(fn)):
			if fn[-(i+1)]=='/':
				return fn[-i:]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def extract_path(self,fn):
		for i in range(len(fn)):
			if fn[-(i+1)]=='/':
				return fn[:-i]
#                              MISC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def ExtractRotMat(mat):#from 4x4 vtk pos/ori matrix
	res=[[0,0,0],[0,0,0],[0,0,0]]
	for i in range(3):
		for j in range(3):
			res[i][j]=mat.GetElement(i,j)
	return res
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def almost(a,b,tolerance=1e-7):
 	return ( abs(a-b)<tolerance )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def R2Eul(R,tolerance=1e-5):
	"""
	R must be an indexable of shape (3,3) and represent and ORTHOGONAL POSITIVE
	DEFINITE matrix. Otherwise, one should use GetEul().
	"""
	fuzz=1e-3
	R=asarray(R,float)
	if det(R) < 0. :
		raise Exception, "determinant is negative\n"+str(R)
	if not allclose(mat(R)*R.T,identity(3),atol=tolerance):
		raise Exception, "not an orthogonal matrix\n"+str(R)
	cang = 2.0-sum(square([ R[0,2],R[1,2],R[2,0],R[2,1],R[2,2] ]))
	cang = sqrt(min(max(cang,0.0),1.0))
	if (R[2,2]<0.0): cang=-cang
	ang=acos(cang)
	beta=degrees(ang)
	sang=sin(ang)
	if(sang>fuzz):
		alpha=degrees(atan2(R[1,2], R[0,2]))
		gamma=degrees(atan2(R[2,1],-R[2,0]))
	else:
		alpha=degrees(atan2(-R[0,1],R[0,0]*R[2,2]))
		gamma=0.
	if   almost(beta,0.,fuzz):
		alpha,beta,gamma = alpha+gamma,  0.,0.
	elif almost(beta,180.,fuzz):
		alpha,beta,gamma = alpha-gamma,180.,0.
	alpha=mod(alpha,360.);
	gamma=mod(gamma,360.)
	if almost(alpha,360.,fuzz):
		alpha=0.
	if almost(gamma,360.,fuzz):
		gamma=0.
	return alpha,beta,gamma
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def Eul2R(alpha,beta,gamma):
  	"""
	Returns the rotation matrix given its eulerian representation (ZYZ
	convention). Angles in degrees.
	"""
	a,b,g = radians(alpha), radians(beta), radians(gamma)
	Ra,Rb,Rg        = [ mat(identity(3))  for i in range(3) ]
	Ra[0:2,0:2]     = [[cos(a),-sin(a)],[sin(a),cos(a)]]
	Rg[0:2,0:2]     = [[cos(g),-sin(g)],[sin(g),cos(g)]]
	Rb[0:3:2,0:3:2] = [[cos(b),sin(b)],[-sin(b),cos(b)]]
	return Ra*Rb*Rg
#===============================================================================
class ps(vtk.vtkPlatonicSolidSource):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def __init__(self,gfx,solidtype,radius,ori,cTU,delta,henttype,dnn,cnn,axe,symfilevar,phizero,caller):
		self.gfx = gfx
		if radius != None:
			if radius.get().rstrip()=='':
				MB.showwarning('Info','Enter a value for the Support Radius')
				return
			else :
				self.radius = float(radius.get())
				if gfx.map!=[]:
					self.radius = self.radius/gfx.map[0].scale
		else:
			self.radius = None
		self.acteur=None
		self.sphact= None
		self.enttype=None
		self.mapper=None
		self.spheremapper=None
		self.solidtype=solidtype
		self.henttype=henttype #helical entry type : ctu or helm-hel
		self.deltas = delta
		if phizero != None:
			self.phizero=float(phizero)
		else :
			self.phizero=0
		self.c=''
		self.T=''
		self.U=''
		self.dphi=''
		self.dz=''
		self.s=1
		self.dnn=1
		self.cnn=1
		self.axe='Z'
		self.ori=ori
		self.numpts=None
		self.symfile=None
		if solidtype=='Icosahedral':
			if self.ori != '':
				self.SetSolidTypeToIcosahedron()
				self.display_platonic(gfx,ori)
			else :
				MB.showwarning('Info','Configure symmetry setting')
		elif solidtype=='Tetrahedral':
			if self.ori != '':
				self.SetSolidTypeToTetrahedron()
				self.display_platonic(gfx,ori)
			else :
				MB.showwarning('Info','Configure symmetry setting')
		elif solidtype=='Octahedral':
			if self.ori != '':
				self.SetSolidTypeToOctahedron()
				self.display_platonic(gfx,ori)
			else :
				MB.showwarning('Info','Configure symmetry setting')
		elif solidtype=='Helicoidal':
			if cTU[0].get()!='':
				self.c=float(cTU[0].get())
			if cTU[1].get()!='':
				self.T=int(float(cTU[1].get()))
			if cTU[2].get()!='':
				self.U=int(float(cTU[2].get()))
			if delta[0].get()!='':
				self.dphi=float(delta[0].get())
			if delta[1].get()!='':
				self.dz=float(delta[1].get())
			if delta[2].get()!='':
				self.s=int(delta[2].get())
			if not((self.c!='' and self.T!='' and self.U!='') or (self.dphi!='' and self.dz!='')):
				MB.showwarning('Info','Configure symmetry setting')
				return
			self.display_tube(gfx,caller='sym')
			return
		elif solidtype=='Dn':
			if dnn != '' :
				if gfx.ps != None:
					gfx.renderer.RemoveActor(gfx.ps.acteur)
					gfx.renwin.Render()
				self.dnn = int(dnn)
				system(self.gfx.vedabin + '/symmetry.exe <<ENDOF\n dn \n %d \nENDOF'%self.dnn)
				self.display_Xn(gfx)
				print 'the Symmetry D%d has been defined'%self.dnn
				gfx.ps = self
		elif solidtype=='Cn':
			if cnn != '' :
				if gfx.ps != None:
					gfx.renderer.RemoveActor(gfx.ps.acteur)
					gfx.renwin.Render()
				self.cnn = int(cnn)
				self.axe=axe
				system(self.gfx.vedabin + '/symmetry.exe <<ENDOF\n cn \n %d %s \nENDOF'%(self.cnn,self.axe))
				self.display_Xn(gfx)
				print 'the Symmetry C%d has been defined around axe %s'%(self.cnn,self.axe)
				gfx.ps = self
		elif solidtype=='P1':
			if gfx.ps != None:
				gfx.renderer.RemoveActor(gfx.ps.acteur)
				gfx.renwin.Render()
			system(self.gfx.vedabin + '/symmetry.exe >> /dev/null <<ENDOF\n p1 \nENDOF')
			print 'the Symmetry P1 has been defined'
			gfx.ps = self
		elif solidtype=='User Defined':
			if symfilevar.get() != '':
				if gfx.ps != None:
					gfx.renderer.RemoveActor(gfx.ps.acteur)
				self.symfile = symfilevar.get()
				system('cp %s %s/SYM'%(self.symfile,self.gfx.workdir))
				print 'the "User Defined" symmetry has been defined'
				gfx.ps = self
		else :
			print 'Option not yet implemented'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def display_Xn(self,gfx):
		if gfx.ps != None:
			gfx.renderer.RemoveActor(gfx.ps.acteur)
		try:
			zmax = gfx.map[0].box.GetBounds()[5]
			zmin = gfx.map[0].box.GetBounds()[4]
			ymax = gfx.map[0].box.GetBounds()[3]
			ymin = gfx.map[0].box.GetBounds()[2]
			c = (zmax - zmin)/2.
			b = (ymax - ymin)/2.
		except:
				b = c = self.radius
		dist = self.radius
		phizero = radians(self.phizero)
		if self.solidtype=='Cn':
			numPts = self.cnn
			shifta = radians(360./numPts)
			pdo = vtk.vtkPolyData()
			cirpdo = vtk.vtkPolyData()
			newPts = vtk.vtkPoints()
			vals = vtk.vtkDoubleArray()
			for i in range(0, numPts):
				try :
					if self.axe=='Z':
						x = dist * gfx.map[0].scale * cos(i*shifta + phizero)
						y = dist * gfx.map[0].scale * sin(i*shifta + phizero)
						z = 0
					elif self.axe == 'Y':
						x = dist * gfx.map[0].scale * cos(i*shifta + phizero)
						y = 0
						z = dist * gfx.map[0].scale * sin(i*shifta + phizero)
					elif self.axe == 'X':
						x = 0
						y = dist * gfx.map[0].scale * cos(i*shifta + phizero)
						z = dist * gfx.map[0].scale * sin(i*shifta + phizero)
				except :
					if self.axe=='Z':
						x = dist  * cos(i*shifta + phizero)
						y = dist  * sin(i*shifta + phizero)
						z = 0
					elif self.axe == 'Y':
						x = dist  * cos(i*shifta + phizero)
						y = 0
						z = dist  * sin(i*shifta + phizero)
					elif self.axe == 'X':
						x = 0
						y = dist * cos(i*shifta + phizero)
						z = dist * sin(i*shifta + phizero)
				newPts.InsertPoint(i, x,y,z)
				vals.InsertNextValue(i)
			pdo.SetPoints(newPts)
			pdo.GetPointData().SetScalars(vals)
			cirpdo.SetPoints(newPts)
			cirpdo.GetPointData().SetScalars(vals)
			aPolyLine = vtk.vtkPolyLine()
			aPolyLine.GetPointIds().SetNumberOfIds(numPts+1)
			for i in range(0,numPts):
				aPolyLine.GetPointIds().SetId(i, i)
			aPolyLine.GetPointIds().SetId(numPts, 0)
			cirpdo.Allocate(1, 1)
			cirpdo.InsertNextCell(aPolyLine.GetCellType(), aPolyLine.GetPointIds())
			sphere = vtk.vtkSphereSource()
			sphere.SetCenter(0, 0, 0)
			try :
				sphere.SetRadius(self.radius*gfx.map[0].scale/5.0)
			except:
				sphere.SetRadius(self.radius/5.0)
			sphere.SetThetaResolution(8)
			sphere.SetStartTheta(0)
			sphere.SetEndTheta(360)
			sphere.SetPhiResolution(8)
			sphere.SetStartPhi(0)
			sphere.SetEndPhi(180)
			glyph = vtk.vtkGlyph3D()
			glyph.SetInput(pdo)
			glyph.SetColorMode(1)
			glyph.ScalingOn()
			glyph.SetScaleMode(2)
			glyph.SetScaleFactor(0.25)
			glyph.SetSource(sphere.GetOutput())
			self.spheremapper = vtk.vtkPolyDataMapper()
			self.spheremapper.SetInputConnection(glyph.GetOutputPort())
			self.spheremapper.UseLookupTableScalarRangeOff()
			self.spheremapper.SetScalarVisibility(0)
			self.sphact = vtk.vtkActor()
			self.sphact.PickableOff()
			self.sphact.DragableOff()
			self.sphact.SetMapper(self.spheremapper)
			self.sphact.GetProperty().SetRepresentationToSurface()
			self.sphact.GetProperty().SetInterpolationToGouraud()
			self.sphact.GetProperty().SetAmbient(0.15)
			self.sphact.GetProperty().SetDiffuse(0.85)
			self.sphact.GetProperty().SetSpecular(0.1)
			self.sphact.GetProperty().SetSpecularPower(100)
			self.sphact.GetProperty().SetSpecularColor(1, 1, 1)
			self.sphact.GetProperty().SetColor(1,1,1)
			self.mapper = vtk.vtkPolyDataMapper()
			self.mapper.SetScalarVisibility(0)
			self.mapper.SetInput(cirpdo)
			self.helact = vtk.vtkActor()
			self.helact.PickableOff()
			self.helact.DragableOff()
			self.helact.GetProperty().SetColor(1,1,1)
			self.helact.GetProperty().SetRepresentationToWireframe()
			self.helact.GetProperty().SetLineWidth(5)
			self.helact.GetProperty().SetLineWidth(1)
			self.helact.GetProperty().SetSpecular(.4)
			self.helact.GetProperty().SetSpecularPower(10)
			self.helact.SetMapper(self.mapper)
			assembly = vtk.vtkAssembly()
			assembly.AddPart(self.sphact)
			assembly.AddPart(self.helact)
			self.acteur=assembly
			gfx.renderer.AddActor(self.acteur)
			gfx.renwin.Render()
			gfx.ps = self
		elif self.solidtype=='Dn':
			numPts = self.dnn
			shifta = radians(360./numPts)
			cirpdo = vtk.vtkPolyData()
			pdo = vtk.vtkPolyData()
			newPts = vtk.vtkPoints()
			vals = vtk.vtkDoubleArray()
			for i in range(0, numPts):
				try:
					x = dist * gfx.map[0].scale * cos(i*shifta)
					y = dist * gfx.map[0].scale * sin(i*shifta)
					z = 0
				except:
					x = dist * cos(i*shifta)
					y = dist * sin(i*shifta)
					z = 0
				newPts.InsertPoint(i, x,y,z)
				vals.InsertNextValue(i)
			pdo.SetPoints(newPts)
			pdo.GetPointData().SetScalars(vals)
			cirpdo.SetPoints(newPts)
			cirpdo.GetPointData().SetScalars(vals)
			aPolyLine = vtk.vtkPolyLine()
			aPolyLine.GetPointIds().SetNumberOfIds(numPts+1)
			for i in range(0,numPts):
				aPolyLine.GetPointIds().SetId(i, i)
			aPolyLine.GetPointIds().SetId(numPts, 0)
			cirpdo.Allocate(1, 1)
			cirpdo.InsertNextCell(aPolyLine.GetCellType(), aPolyLine.GetPointIds())
			sphere = vtk.vtkSphereSource()
			sphere.SetCenter(0, 0, 0)
			try:
				sphere.SetRadius(self.radius*gfx.map[0].scale/5.0)
			except:
				sphere.SetRadius(self.radius/5.0)
			sphere.SetThetaResolution(8)
			sphere.SetStartTheta(0)
			sphere.SetEndTheta(360)
			sphere.SetPhiResolution(8)
			sphere.SetStartPhi(0)
			sphere.SetEndPhi(180)
			glyph = vtk.vtkGlyph3D()
			glyph.SetInput(pdo)
			glyph.SetColorMode(1)
			glyph.ScalingOn()
			glyph.SetScaleMode(2)
			glyph.SetScaleFactor(0.25)
			glyph.SetSource(sphere.GetOutput())
			self.spheremapper = vtk.vtkPolyDataMapper()
			self.spheremapper.SetInputConnection(glyph.GetOutputPort())
			self.spheremapper.UseLookupTableScalarRangeOff()
			self.spheremapper.SetScalarVisibility(0)
			self.sphact = vtk.vtkActor()
			self.sphact.PickableOff()
			self.sphact.DragableOff()
			self.sphact.SetMapper(self.spheremapper)
			self.sphact.GetProperty().SetRepresentationToSurface()
			self.sphact.GetProperty().SetInterpolationToGouraud()
			self.sphact.GetProperty().SetAmbient(0.15)
			self.sphact.GetProperty().SetDiffuse(0.85)
			self.sphact.GetProperty().SetSpecular(0.1)
			self.sphact.GetProperty().SetSpecularPower(100)
			self.sphact.GetProperty().SetSpecularColor(1, 1, 1)
			self.sphact.GetProperty().SetColor(1,1,1)
			self.mapper = vtk.vtkPolyDataMapper()
			self.mapper.SetScalarVisibility(0)
			self.mapper.SetInput(cirpdo)
			self.helact = vtk.vtkActor()
			self.helact.PickableOff()
			self.helact.DragableOff()
			self.helact.GetProperty().SetColor(1,1,1)
			self.helact.GetProperty().SetRepresentationToWireframe()
			self.helact.GetProperty().SetLineWidth(5)
			self.helact.GetProperty().SetLineWidth(1)
			self.helact.GetProperty().SetSpecular(.4)
			self.helact.GetProperty().SetSpecularPower(10)
			self.helact.SetMapper(self.mapper)
			sphact2 = vtk.vtkActor()
			helact2 = vtk.vtkActor()
			sphact2.ShallowCopy(self.sphact)
			helact2.ShallowCopy(self.helact)
			self.sphact.AddPosition(0,0,dist*c/b)
			self.helact.AddPosition(0,0,dist*c/b)
			sphact2.AddPosition(0,0,-dist*c/b)
			helact2.AddPosition(0,0,-dist*c/b)
			assembly = vtk.vtkAssembly()
			assembly.AddPart(self.sphact)
			assembly.AddPart(self.helact)
			assembly.AddPart(sphact2)
			assembly.AddPart(helact2)
			self.acteur=assembly
			gfx.renderer.AddActor(self.acteur)
			gfx.renwin.Render()
			gfx.ps = self
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def display_tube(self,gfx,caller):
		if gfx.ps != None:
			gfx.renderer.RemoveActor(gfx.ps.acteur)
		if self.c!='' and self.T!='' and self.U!='' and self.henttype=='cTU':#CTU CASE
			ptype='delta'
			self.enttype='cTU'
			deltaz = (self.c/self.U)
			deltaphi = (self.T * 360.0)/self.U
			self.dphi = deltaphi
			self.dz = deltaz
		elif self.dphi!='' and self.dz!='' and self.henttype=='Elm-Hel':#DELTA CASE
			ptype = 'delta'
			self.enttype='Elm-Hel'
			deltaphi = self.dphi
			deltaz =self.dz
		else :
			MB.showwarning('Info','Configure symmetry setting')
			return
		try:
			zmax = gfx.map[0].box.GetBounds()[5]
			zmin = gfx.map[0].box.GetBounds()[4]
			z = zmax - zmin
		except:
			z= 4*360*abs(self.dz/self.dphi)
		#ici exeption 'fit' calcule de dz avec ratio !!!
		if 'fit' in caller:
			self.c = self.c * gfx.map[0].ratio
			self.dz = self.dz * gfx.map[0].ratio
			deltaz =self.dz
		numPts = int(z/deltaz)
		self.numpts = numPts#definition du nombre de sm
		self.cptsymops_tube(ptype,numPts) #ptype est tjs = a delta on calcule toujours de la même maniere dans symmetry
		#parameter definition :
		shifta = radians(deltaphi)
		phizero = radians(self.phizero)
		###definition des boules selon les parametres
		pdo = vtk.vtkPolyData()
		newPts = vtk.vtkPoints() #This will store the points for the Helix
		vals = vtk.vtkDoubleArray()
		for i in range(0, numPts):
			try :
				x = self.radius * gfx.map[0].scale * cos(i*shifta + phizero)
				y = self.radius * gfx.map[0].scale * sin(i*shifta + phizero)
				z = zmin + i*deltaz
			except :
				x = self.radius *  cos(i*shifta + phizero)
				y = self.radius * sin(i*shifta + phizero)
				z =  i*deltaz
			newPts.InsertPoint(i, x,y,z)
			vals.InsertNextValue(i)
		pdo.SetPoints(newPts)
		pdo.GetPointData().SetScalars(vals)
		sphere = vtk.vtkSphereSource()
		sphere.SetCenter(0, 0, 0)
		try :
			sphere.SetRadius(self.radius*gfx.map[0].scale/5.0)
		except:
			sphere.SetRadius(self.radius*(1/5.))
		sphere.SetThetaResolution(8)
		sphere.SetStartTheta(0)
		sphere.SetEndTheta(360)
		sphere.SetPhiResolution(8)
		sphere.SetStartPhi(0)
		sphere.SetEndPhi(180)
		glyph = vtk.vtkGlyph3D()
		glyph.SetInput(pdo)
		glyph.SetColorMode(1)
		glyph.ScalingOn()
		glyph.SetScaleMode(2)
		glyph.SetScaleFactor(0.25)
		glyph.SetSource(sphere.GetOutput())
		self.spheremapper = vtk.vtkPolyDataMapper()
		self.spheremapper.SetInputConnection(glyph.GetOutputPort())
		self.spheremapper.UseLookupTableScalarRangeOff()
		self.spheremapper.SetScalarVisibility(0)
		self.sphact = vtk.vtkActor()
		self.sphact.PickableOff()
		self.sphact.DragableOff()
		self.sphact.SetMapper(self.spheremapper)
		self.sphact.GetProperty().SetRepresentationToSurface()
		self.sphact.GetProperty().SetInterpolationToGouraud()
		self.sphact.GetProperty().SetAmbient(0.15)
		self.sphact.GetProperty().SetDiffuse(0.85)
		self.sphact.GetProperty().SetSpecular(0.1)
		self.sphact.GetProperty().SetSpecularPower(100)
		self.sphact.GetProperty().SetSpecularColor(1, 1, 1)
		self.sphact.GetProperty().SetColor(1,1,1)
		###definition de l helice continue
		enhance=100
		hlxpdo = vtk.vtkPolyData()
		hlxnewPts = vtk.vtkPoints() #This will store the points for the Helix
		hlxvals = vtk.vtkDoubleArray()
		hlxnumPts = int(numPts*enhance)
		for i in range(0, hlxnumPts):
			try :
				x = self.radius * gfx.map[0].scale * cos(i*shifta/enhance + phizero)
				y = self.radius * gfx.map[0].scale * sin(i*shifta/enhance + phizero)
				z = zmin + i*deltaz/enhance
			except :
				x = self.radius  * cos(i*shifta/enhance + phizero)
				y = self.radius  * sin(i*shifta/enhance + phizero)
				z = i*deltaz/enhance
			hlxnewPts.InsertPoint(i, x,y,z)
			hlxvals.InsertNextValue(i)
		hlxpdo.SetPoints(hlxnewPts)
		hlxpdo.GetPointData().SetScalars(hlxvals)
		aPolyLine = vtk.vtkPolyLine()
		aPolyLine.GetPointIds().SetNumberOfIds(hlxnumPts)
		for i in range(0,hlxnumPts):
			aPolyLine.GetPointIds().SetId(i, i)
		hlxpdo.Allocate(1, 1)
		hlxpdo.InsertNextCell(aPolyLine.GetCellType(), aPolyLine.GetPointIds())
		self.mapper = vtk.vtkPolyDataMapper()
		self.mapper.SetScalarVisibility(0)
		self.mapper.SetInput(hlxpdo)
		self.helact = vtk.vtkActor()
		self.helact.PickableOff()
		self.helact.DragableOff()
		self.helact.GetProperty().SetColor(1,1,1)
		self.helact.GetProperty().SetRepresentationToWireframe()
		self.helact.GetProperty().SetLineWidth(5)
		self.helact.GetProperty().SetLineWidth(1)
		self.helact.GetProperty().SetSpecular(.4)
		self.helact.GetProperty().SetSpecularPower(10)
		self.helact.SetMapper(self.mapper)
		assembly = vtk.vtkAssembly()
		assembly.AddPart(self.sphact)
		assembly.AddPart(self.helact)
		helcol = [(1,1,1),(1,0,0.3),(0.3,0,1),(0,1,0.3)]
		for i in range(1,self.s+1):
			if i == 1:
				continue
			else :
				tmp_sphact= vtk.vtkActor()
				tmp_helact= vtk.vtkActor()
				tmp_sphact.ShallowCopy(self.sphact)
				tmp_helact.ShallowCopy(self.helact)
				ps=vtk.vtkProperty()
				ph=vtk.vtkProperty()
				ps.SetColor(helcol[(i-1)%4])
				ph.SetColor(helcol[(i-1)%4])
				tmp_sphact.SetProperty(ps)
				tmp_helact.SetProperty(ph)
				tmp_sphact.RotateWXYZ((360/self.s)*(i-1),0,0,1)
				tmp_helact.RotateWXYZ((360/self.s)*(i-1),0,0,1)
				assembly.AddPart(tmp_sphact)
				assembly.AddPart(tmp_helact)
		self.acteur=assembly
		gfx.renderer.AddActor(self.acteur)
		gfx.renwin.Render()
		gfx.ps = self
		return ptype
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def cptsymops_tube(self,ptype,numPts):
		if ptype == 'ctu':#n est pas utilise
			system(self.gfx.vedabin + '/symmetry.exe <<ENDOF\n hel \n ctu \n %f %d %d %d \nENDOF'%(self.c,self.T,self.U,numPts))
		elif ptype == 'delta':
			system(self.gfx.vedabin + '/symmetry.exe <<ENDOF\n hel \n ele \n %f %f %d %d  \nENDOF'%(self.dphi,self.dz,self.s,1.5*numPts))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def display_platonic(self,gfx,ori):
		if gfx.ps != None:
			gfx.renderer.RemoveActor(gfx.ps.acteur)
		self.mapper = vtk.vtkPolyDataMapper()
		self.mapper.SetInputConnection(self.GetOutputPort())
		self.mapper.SetScalarVisibility(0)
		self.acteur = vtk.vtkActor()
		self.acteur.PickableOff()
		self.acteur.SetMapper(self.mapper)
		if self.solidtype=='Icosahedral':
			self.orient_ico(ori)
		elif self.solidtype=='Octahedral':
			self.orient_octa(ori)
		elif self.solidtype=='Tetrahedral':
			self.orient_tetra(ori)
		self.acteur.GetProperty().SetColor(1,1,1)
		try:
			self.acteur.SetScale(self.radius*gfx.map[0].scale,self.radius*gfx.map[0].scale,self.radius*gfx.map[0].scale)
		except :
			self.acteur.SetScale(self.radius,self.radius,self.radius)
		self.acteur.GetProperty().SetRepresentationToWireframe()
		self.acteur.GetProperty().SetInterpolationToPhong()
		self.acteur.GetProperty().SetLineWidth(5)
		self.acteur.GetProperty().ShadingOn()
		self.acteur.GetProperty().SetAmbient(0.6)
		self.acteur.GetProperty().SetDiffuse(0.4)
		gfx.renderer.AddActor(self.acteur)
		gfx.renwin.Render()
		gfx.ps = self
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def orient_ico(self,ori):
		pid=degrees(pi)
		pio2d=degrees(pi/2.)
		ang2=degrees(acos(1./sqrt(5.)))/2.
		ang3=degrees(acos(sin(pi*3./5.)*(1.+sqrt(5.))/sqrt(15.))) #il faut passer a radians
		zero=0.
		#from 2Z2Y.2 to 5Z2YNX :
		self.RotaEuler(zero,pio2d-ang2,zero)
		if   ori=='5Z2Y.1':
			pass
		elif ori=='5Z2Y.2':
			self.RotaEuler(pid,pid,zero)
		elif ori=='5Z2X.1':
			self.RotaEuler(pio2d,zero,zero)
		elif ori=='5Z2X.2':
			self.RotaEuler(-pio2d,zero,zero)
		elif ori=='5Y2Z.1':
			self.RotaEuler(pio2d,pio2d,pio2d)
		elif ori=='5Y2Z.2':
			self.RotaEuler(pio2d,pio2d,-pio2d)
		elif ori=='5Y2X.1':
			self.RotaEuler(pio2d,pio2d,pid)
		elif ori=='5Y2X.2':
			self.RotaEuler(pio2d,pio2d,zero)
		elif ori=='5X2Z.1':
			self.RotaEuler(zero,pio2d,pio2d)
		elif ori=='5X2Z.2':
			self.RotaEuler(zero,pio2d,-pio2d)
		elif ori=='5X2Y.1':
			self.RotaEuler(zero,pio2d,pid)
		elif ori=='5X2Y.2':
			self.RotaEuler(zero,pio2d,zero)
		elif ori=='3Z2Y.1':
			self.RotaEuler(pid,ang3,pid)
		elif ori=='3Z2Y.2':
			self.RotaEuler(zero,ang3,pid)
		elif ori=='3Z2X.1':
			self.RotaEuler(pio2d,ang3,pid)
		elif ori=='3Z2X.2':
			self.RotaEuler(-pio2d,ang3,pid)
		elif ori=='3Y2Z.1':
			self.RotaEuler(pio2d-ang3,pio2d,pio2d)
		elif ori=='3Y2Z.2':
			self.RotaEuler(3*pio2d+ang3,pio2d,-pio2d)
		elif ori=='3Y2X.1':
			self.RotaEuler(pio2d,pio2d+ang3,pid)
		elif ori=='3Y2X.2':
			self.RotaEuler(pio2d,pio2d-ang3,zero)
		elif ori=='3X2Z.1':
			self.RotaEuler(pi-ang3,pio2d,pio2d)
		elif ori=='3X2Z.2':
			self.RotaEuler(ang3,pio2d,-pio2d)
		elif ori=='3X2Y.1':
			self.RotaEuler(pid,pio2d-ang3,zero)
		elif ori=='3X2Y.2':
			self.RotaEuler(zero,pio2d-ang3,zero)
		elif ori=='2Z2Y.1':
			self.RotaEuler(zero,ang2,zero)
		elif ori=='2Z2Y.2':
			self.RotaEuler(pid,pio2d-ang2,pid)
		else : print 'another sym option'
		system(self.gfx.vedabin + '/symmetry.exe <<ENDOF\n ico \n%s \nENDOF'%ori)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def orient_octa(self,ori):
		pid=degrees(pi)
		pio2d=degrees(pi/2.)
		pio3d=degrees(pi/3.)
		pio4d=degrees(pi/4.)
		ang2=degrees(acos(1./sqrt(3.)))
		zero=0.
		#from 4Z2Y to 4Z4Y4X :
		self.RotaEuler(zero,zero,pio4d)
		if   ori== '4Z4Y':
			pass
		elif ori== '4Z2Y':
			self.RotaEuler(zero,zero,pio4d)
		elif ori== '4Y2Z':
			self.RotaEuler(zero,pio4d,zero)
		elif ori== '4X2Z':
			self.RotaEuler(-pio2d,pio4d,zero)
		elif ori== '3Z2Y.1':
			self.RotaEuler(zero,ang2,3*pio4d)
		elif ori== '3Z2Y.2':
			self.RotaEuler(pio3d,ang2,3*pio4d)
		elif ori== '3Z2X.1':
			self.RotaEuler(pio2d,ang2,3*pio4d)
		elif ori== '3Z2X.2':
			self.RotaEuler(-pio2d,ang2,3*pio4d)
		elif ori== '3Y2Z.1':
			self.RotaEuler(pio2d-ang2,pio2d,pio4d)
		elif ori== '3Y2Z.2':
			self.RotaEuler(ang2,3*pio4d,pio2d)
		elif ori== '3Y2X.1':
			self.RotaEuler(pio2d,pio2d-ang2,-pio4d)
		elif ori== '3Y2X.2':
			self.RotaEuler(pio2d,pio2d+ang2,3*pio4d)
		elif ori== '3X2Z.1':
			self.RotaEuler(ang2,pio2d,5*pio4d)
		elif ori== '3X2Z.2':
			self.RotaEuler(pio2d-ang2,pio4d,-pio2d)
		elif ori== '3X2Y.1':
			self.RotaEuler(zero,pio2d-ang2,-pio4d)
		elif ori== '3X2Y.2':
			self.RotaEuler(zero,pio2d+ang2,3*pio4d)
		else : print 'another sym option'
		system(self.gfx.vedabin + '/symmetry.exe <<ENDOF\n oct \n%s \nENDOF'%ori)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def orient_tetra(self,ori):
		#REF 2Z2Y2X :
		if   ori==  '2Z2Y2X':
			pass
		system(self.gfx.vedabin + '/symmetry.exe <<ENDOF\n tet \n%s \nENDOF'%ori)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def RotaEuler(self,alpha,beta,gamma):
		self.acteur.RotateWXYZ(gamma,0,0,1)
		self.acteur.RotateWXYZ(beta,0,1,0)
		self.acteur.RotateWXYZ(alpha,0,0,1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def browsesymfile(gfx,symfilevar):
	fn = FD.askopenfilename(title="Open Sym file",filetypes=[("All", "*")])
	symfilevar.set(fn)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def hide_show_symsupport(gfx):
	if gfx.ps:
		if gfx.ps.acteur.GetVisibility():
			gfx.ps.acteur.VisibilityOff()
		else:
			gfx.ps.acteur.VisibilityOn()
		gfx.renwin.Render()
	else :
		MB.showwarning('Info','Symmetry Support does not Exist')
