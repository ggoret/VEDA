# -*- coding: utf-8 -*-
# URO Fitting Module of VEDA
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
from subprocess import *
from numpy import *
import gfx ,Map,mod, itf,sym
from os import chdir,system, environ,path,putenv
import tkMessageBox as MB
#===============================================================================
class Ifitting():
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def __init__(self,gfx,res):
		self.p=None
		self.itfo=gfx.itf
		self.gfx=gfx
		self.res = res
		gfx.ifit = self
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def Ifitin_init(self):
		self.itfo.ccbar.show()
		self.p = Popen(['e/fitin', 'G0'],bufsize=0,stdin=PIPE,stdout=PIPE,stderr=PIPE,universal_newlines=True)
		self.change_mol()
		print 'Interactive Fitting is Working'
		self.res.set('Interactive Fitting is Working')
		self.itfo.fitreshighentry.configure(state=tk.DISABLED)
		self.itfo.fitreslowentry.configure(state=tk.DISABLED)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def fitin_io(self,obj=None,eve=None):
		frame = self.frame()
		try :
			self.p.stdin.write(frame)
			l = self.p.stdout.read(77)
			cc = float(l.split()[-1])
			self.itfo.ccbar.update(cc)
			self.res.set('Correlation : %s'%cc)
		except :
			MB.showwarning('Error','Problem in function ifit.fitin_io, process fitin must be dead ... RIP')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def change_mol(self):
		try :
			self.p.stdin.write('#%2d \n'%self.gfx.whichmol)
			l = self.p.stdout.read(31)
		except :
			MB.showwarning('Error','Problem in function ifit.change_mol, process fitin must be dead ... RIP')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def terminate_process(self):
		try :
			self.p.stdin.write('\n')
		except :
			MB.showwarning('Error','Interactive Fitting terminated badly')
			pass
		self.itfo.ccbar.hide()
		self.itfo.ccbar.update('CC')
		self.gfx.ifit=None
		self.itfo.fitreshighentry.configure(state=tk.NORMAL)
		self.itfo.fitreslowentry.configure(state=tk.NORMAL)
		self.res.set('Interactive Fitting Stopped')
		print 'Interactive Fitting Stopped'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def frame(self):
		mol = self.gfx.mol[self.gfx.whichmol-1]
		mat=mol.acteur.GetMatrix()
		(a,b,g) = sym.R2Eul(sym.ExtractRotMat(mat))
		(x,y,z) = [mat.GetElement(0,3),mat.GetElement(1,3),mat.GetElement(2,3)]
		return '%.3f %.3f %.3f %.3f %.3f %.3f \n'%(a,b,g,x,y,z)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def interactive_fitting(gfx,reslow,reshigh,res):
	if gfx.mol!=[]:
		if gfx.fit!=None and fit_is_possible(gfx):
			if sym_is_ok(gfx):
				if gfx.ifit!=None:
					gfx.ifit.terminate_process()
				else :
					init_ifit_file(gfx,reslow,reshigh)
					ifit=Ifitting(gfx,res)
					ifit.Ifitin_init()
			else :
				MB.showwarning('Info','There is no sym-mate molecule within the volume. Update Sym-Mates')
		else :
			MB.showwarning('Info','Configure Fitting|Setup')
	else :
		MB.showwarning('Info','Go to Assemblage|Molecules and load molecules. Then configure Fitting|Setup')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def init_ifit_file(gfx,reslow,reshigh):
	gfx.fit.build_fitin_file()
	reso = open('i/reso.i','w')
	reso.write(' %s %s\n'%(reslow,reshigh))
	reso.close()
	system('e/udi2i')
	for mod in gfx.mod:
		if mod.fn in gfx.fit.modlist.keys():
			system('e/tab2i %d'%mod.id)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def reset_ifit(gfx):
	res = gfx.ifit.res
	gfx.ifit.terminate_process()
	gfx.fit.build_fitin_file()
	ifit=Ifitting(gfx,res)
	ifit.Ifitin_init()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def ifit_restart_needed(gfx):
	if gfx.ifit != None:
		gfx.ifit.terminate_process()
		MB.showwarning('Warning','Interactive Fitting has been stopped')
#===============================================================================
#                              URO Refinement fitting
#===============================================================================
class fitting():
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def __init__(self,gfx,reslow,reshigh,vartest):
		self.gfx=gfx
		self.reslow=float(reslow)
		self.reshigh=float(reshigh)
		self.vartest=int(vartest)
		self.fitreslow=float(reslow)
		self.fitreshigh=float(reshigh)
		self.nbc=len(gfx.mol)
		self.nbit=10
		self.varbfact=0
		self.nbmodinuse = 0
		self.modlist=None
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def setup_uro_refinement(self):
		if self.gfx.map==[]:
			MB.showwarning('Info','Load a map')
			self.gfx.itf.status.clear()
			self.gfx.itf.root.configure(cursor='arrow')
			return
		if self.gfx.mol==[]:
			MB.showwarning('Info','Load the molecules')
			self.gfx.itf.status.clear()
			self.gfx.itf.root.configure(cursor='arrow')
			return
		if self.gfx.ps==None:
			MB.showwarning('Info','Define the symmetry (also P1)')
			self.gfx.itf.status.clear()
			self.gfx.itf.root.configure(cursor='arrow')
			return
		if self.gfx.ps.solidtype=='P1':
			self.gfx.sym=sym.Sym(None,self.gfx,0.,0,0)
			self.gfx.sym.refresh_all_constelation(self.gfx)
			self.gfx.itf.root.configure(cursor='watch')
			self.gfx.itf.status.set('Configuration in progress ... please wait')
		if self.vartest==1:
			test = 'test'
		else :
			test = ''
		chdir(self.gfx.workdir)
		param = self.seturo_command_line()
		if self.gfx.sym==None:
			MB.showwarning('Info','Go to Assemblage|Molecules|Constellation and build a constellation')
			return
		system('$VEDA/uro/seturo %s %s'%(param,test))
		self.check_setup()
		self.gfx.fit=self
		#self.draw_sfd_hist()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def check_setup(self):
		dim_error_msg = 'Problem during URO setup; go to help'
		if not path.exists('d/data.d'):
			MB.showwarning('Warning',dim_error_msg)
			return
		if not path.exists('o/sort.s'):
			MB.showwarning('Warning',dim_error_msg)
			return
		else :
			try:
				f=open('o/sort.s',"r")
			except IOError:
				MB.showwarning('Warning',dim_error_msg)
			for l in f:
				if 'stop' in l:
					MB.showwarning('Warning',dim_error_msg)
					return
			f.close()
			for mod in self.gfx.mod:
				if mod.fn in self.modlist.keys():
					if not path.exists('o/tabl%d.s'%mod.id):
						MB.showwarning('Warning','model %d '%i + dim_error_msg)
					else :
						f=open('o/tabl%d.s'%mod.id,"r")
						for l in f:
							if 'stop' in l:
								MB.showwarning('Warning',dim_error_msg)
								return
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def seturo_command_line(self):
		"""
		Command ligne builder
		"""
		emap='--emap='+self.gfx.map[0].fn
		model='  --models='
		modlist={}
		for mol in self.gfx.mol:
			if mol.mod.fn in modlist.keys():
				modlist[mol.mod.fn]+=1
			else:
				modlist[mol.mod.fn]=1
		self.nbmodinuse = len(modlist.keys())
		self.modlist = modlist
		for mod in self.gfx.mod:
			if mod.fn in modlist.keys():
				if mod.type=='map':
					model += "%s '#' %d : %d,"%(mod.dfn,mod.id,modlist[mod.fn])
				else :
					model += "%s '#' %d : %d,"%(mod.dfn,mod.id,modlist[mod.fn])
		resol='  --resol=%.2f %.2f'%(self.reslow,self.reshigh)
		scale='  --scale=1.0'#the magnification is now handled in veda itself
		return emap + model + resol + scale
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def draw_sfd_hist(self):
		"""
		Structure factors distribution histogram drawer
		"""
		gp = open('sfd_hist.gp','w')
		gp.write('set title "Wilson plot"\n')
		gp.write('set xtics rotate by 90\n')
		gp.write('set ylabel "ln(<I>)"\n')
		gp.write('set xlabel "Angstroms"\n')
		fmt = "%5.1f -%5.1f"
		col=1
		gp.write("""plot '%s' using 4:xtic(sprintf("%s",$1,$3)) with  linespoints pt 5 lc %d title '%s'\n"""%('HISTU',fmt,col,'Target Map'))
		for mod in self.gfx.mod:
			if mod.fn in self.modlist.keys():
				col+=1
				modf = 'HISTT-%d'%mod.id
				gp.write("""replot '%s' using 4:xtic(sprintf("%s",$1,$3)) with  linespoints pt 5 lc %d title '%s'\n"""%(modf,fmt,col,mod.un))
		gp.close()
		system('gnuplot -persist sfd_hist.gp')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def uro_refinement(self,res):
		if self.fitreshigh < self.reshigh :
			MB.showwarning('Warning','Working resolution range incompatible with Maximal resolution range. Go to Fitting|Setup')
			return
		self.build_fitin_file()
		system('e/docking G0 G1')
		self.assign_pv_from_file('G1',res)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def build_fitin_file(self):
		go = open('G0','w')
		for mol in self.gfx.mol:
			mat=mol.acteur.GetMatrix()
			nbsymop = len(mol.lnbsm)
			strsymlist = str(mol.lnbsm)[1:-1].replace(',','')
			(a,b,g) = sym.R2Eul(sym.ExtractRotMat(mat))
			(x,y,z) = [mat.GetElement(0,3),mat.GetElement(1,3),mat.GetElement(2,3)]
			go.write(' #%2d %.3f %.3f %.3f %.3f %.3f %.3f %d %s\n'%(mol.mod.id,a,b,g,x,y,z,nbsymop,strsymlist))
		go.close()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def assign_pv_from_file(self,g,res):
		try :
			pvf = open('%s'%g,'r')
		except IOError:
			MB.showwarning('Info','Problem during fitting. Go to help')
			return
		molid=self.gfx.mol[0].id
		pvs=[]
		for ligne in pvf:
			if ligne =='':
				MB.showwarning('Info','Problem during fitting. Go to help')
				return
			l=ligne.replace('#','').split()
			if len(l)==9:
				res.set('Correlation : %s, R-factor : %s'%(l[7],l[8]))
			if len(l)>9:
				res.set('Undo mode')
			pvs+=[[molid,float(l[1]),float(l[2]),float(l[3]),float(l[4]),float(l[5]),float(l[6])]]
			molid=self.next_mol_id(molid) #return none if error or end
		for mol in self.gfx.mol:
			for pv in pvs:
				if mol.id==pv[0]:
					self.gfx.sym.render=0 #stop du rendu
					mol.acteur.SetOrientation(0,0,0) #mise a l'origine tres importante
					self.RotaEuler(mol.acteur,pv[1],pv[2],pv[3]) #on defini la partie rotationelle de la transformation
					mol.acteur.SetPosition(pv[4],pv[5],pv[6])
		self.gfx.renwin.Render()
		self.gfx.sym.render=1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def next_mol_id(self,ident):
		for i in range(len(self.gfx.mol)-1):
			if self.gfx.mol[i].id==ident:
				 return	self.gfx.mol[i+1].id
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def build_oicfd(self,fitreslow,fitreshigh,nbc,nbit,varbfact):
		self.fitreslow=fitreslow
		self.fitreshigh=fitreshigh
		self.nbc=nbc
		self.nbit=nbit
		self.varbfact=varbfact
		oic = open('i/oicfd.i2','w')
		oic.write('oic fiting  : Created by VEDA for URO\n')
		oic.write(' ** URO **\n')
		oic.write('  0.0      1000\n')
		oic.write('  %.3f     %.3f RESOLUTION\n'%(fitreslow,fitreshigh))
		oic.write('    %d   %d     0.0    %d CYCLES\n'%(nbc,nbit,varbfact)) #o=condition d'arret
		oic.write('fuzz:       0.0\n')
		oic.write('over:       0.0\n')
		oic.close()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def RotaEuler(self,acteur,alpha,beta,gamma):
		acteur.RotateWXYZ(gamma,0,0,1)
		acteur.RotateWXYZ(beta,0,1,0)
		acteur.RotateWXYZ(alpha,0,0,1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              main functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def configure_fitting(gfx,reslow,reshigh,vartest,scale):
	if reshigh=='' or reslow =='':
		MB.showwarning('Info','Enter a value for Resolution Range')
		return
	if gfx.map==[]:
		MB.showwarning('Info','Load a map')
		return
	gfx.itf.root.configure(cursor='watch')
	gfx.itf.status.set('Configuration in progress ... please wait')
	if scale != gfx.map[0].scale:
			gfx.map[0].rescale_map(gfx,scale,caller='fit')
			gfx.itf.propagate_scale(gfx,scale,caller='fit')
	fit=fitting(gfx,reslow,reshigh,vartest)
	fit.setup_uro_refinement()
	gfx.itf.status.clear()
	gfx.itf.root.configure(cursor='arrow')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def refinement(gfx,fitreslow,fitreshigh,nbc,nbit,varbfact,res):
	if gfx.mol!=[]:
		if gfx.fit!=None and fit_is_possible(gfx):
			if sym_is_ok(gfx):
				gfx.itf.root.configure(cursor='watch')
				gfx.itf.status.set('Refinement in progress ... please wait')
				gfx.fit.build_oicfd(float(fitreslow),float(fitreshigh),int(nbc),int(nbit),int(varbfact))
				gfx.fit.uro_refinement(res)
				gfx.itf.status.clear()
				gfx.itf.root.configure(cursor='arrow')
			else :
				MB.showwarning('Info','There is no sym-mate molecule within the volume. Update Sym-Mates')
		else :
			MB.showwarning('Info','Configure Fitting|Setup')
	else :
		MB.showwarning('Info','Go to Assemblage|Molecules and load molecules. Then configure Fitting|Setup')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def fit_is_possible(gfx):
	chdir(gfx.workdir)
#~~~~~~check map files~~~~~~~~~~~~~~~~~~#
	if not path.exists('d/data.d'):
		return 0
	if not path.exists('f/xudi',):
		return 0
	if not path.exists('o/sort.s'):
		return 0
	else :
		f=open('o/sort.s','r')
		for l in f:
			if 'stop' in l:
				return 0
#~~~~~~check mod files if in uses~~~~~~~#
	modlist={}
	for mol in gfx.mol:
		if mol.mod.fn in modlist.keys():
			modlist[mol.mod.fn]+=1
		else:
			modlist[mol.mod.fn]=1
	for mod in gfx.mod:
		if mod.fn in modlist.keys():
			if not path.exists('o/tabl%d.s'%mod.id) or not path.exists('f/tabl%d'%mod.id):
				return 0
			else :
				f=open('o/tabl%d.s'%mod.id,"r")
				for l in f:
					if 'stop' in l:
						return 0
	return 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def sym_is_ok(gfx):
	for mol in gfx.mol:
		if mol.lnbsm != None:
			if len(mol.lnbsm) == 0:
				return 0
		else :
			return 0
	return 1
#===============================================================================
#                              RMS Threshold
#===============================================================================
class ncs_rms():
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def __init__(self,gfx,rms,stype,nbtrials,nbite):
		self.gfx=gfx
		self.rms=rms
		self.nbtrials=nbtrials
		self.nbite=nbite
		if stype == 'Rotation':
			self.stype='R'
		elif stype == 'Translation':
			self.stype = 'T'
		else :
			self.stype = 'A'
		gfx.ncs_rms=self
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def start(self,mol):
		chdir(self.gfx.workdir)
		self.gfx.fit.build_fitin_file()
		self.gfx.fit.build_oicfd(self.gfx.fit.fitreslow,self.gfx.fit.fitreshigh,self.gfx.fit.nbc,self.nbite,self.gfx.fit.varbfact)
		system('e/ncsrms G0 <<ENDOF & \n%d\n%f\n%s\n%d\nENDOF'%(mol.whoami(self.gfx),self.rms,self.stype,self.nbtrials))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def draw_rms_hist(self):
		"""
		a mettre apres
		"""
		gp = open('rms_hist.gp','w')
		gp.write('set title "RMS convergence threshold"\n')
		gp.write('set xlabel "trial #"\n')
		gp.write("plot 'HISTR' using 10 with linespoints pt 5 lc 1 title 'rms to solution'\n")
		gp.write("replot 'HISTR' using 8 with  linespoints pt 5 lc 2 title '100 x CC'")
		gp.close()
		system('gnuplot -persist rms_hist.gp')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def launch_rms_threshold(gfx,rms,stype,nbtrials,nbite,mollist):
	if gfx.mol!=[]:
		if gfx.fit!=None and fit_is_possible(gfx):
			if sym_is_ok(gfx):
				for mol in gfx.mol:
					if mollist.get(tk.ACTIVE) == mol.un:
						try:
							rms=float(rms)
							nbtrials=int(nbtrials)
							nbite=int(nbite)
						except:
							MB.showwarning('Error','Enter a RMS value')
							return
						rmst=ncs_rms(gfx,rms,stype,nbtrials,nbite)
						rmst.start(mol)
						#rmst.draw_rms_hist()
				gfx.itf.showinfo('Info','URO Refinement is running in background for the different trials')
			else :
				MB.showwarning('Info','There is no sym-mate molecule within the volume. Update Sym-Mates')
		else :
			MB.showwarning('Info','Configure Fitting|Setup')
	else :
		MB.showwarning('Info','Go to Assemblage|Molecules and load molecules. Then configure Fitting|Setup')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def launch_magnum(gfx,step,nbstep):
	if gfx.mol!=[]:
		if gfx.fit!=None and fit_is_possible(gfx):
			if sym_is_ok(gfx):
				step = float(step)
				nbstep = int(nbstep)
				putenv('LC_ALL','C')
				system('e/magnum %f %d &'%(step,nbstep))
				gfx.itf.showinfo('Info','URO Refinement is running in background for the different magnifications')
			else :
				MB.showwarning('Info','There is no sym-mate molecule within the volume. Update Sym-Mates')
		else :
			MB.showwarning('Info','Configure Fitting|Setup')
	else :
		MB.showwarning('Info','Go to Assemblage|Molecules and load molecules. Then configure Fitting|Setup')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_max_rmsd(gfx,mollist):
	chdir(gfx.workdir)
	for mol in gfx.mol:
		try:
			if mollist.get(tk.ACTIVE) == mol.un:
				if path.exists('o/tabl%d.s'%mol.mod.id):
					f=open('o/tabl%d.s'%mol.mod.id,"r")
					for l in f:
						if 'maximal rotational rmsd' in l:
							return l.split()[-1]
					return('NA')
				else :
					return('NA')
		except:
			return('NA')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              CC PROFILER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def cc_profile(gfx,mollist,anglemin,anglemax,step,axe):
	if gfx.ifit!=None:
		MB.showwarning('Info','Stop Interactve Fitting')
		return
	for mol in gfx.mol:
		if mollist.get(tk.ACTIVE) == mol.un:
			chdir(gfx.workdir)
			gfx.fit.build_fitin_file()
			init_ifit_file(gfx,gfx.fit.fitreslow, gfx.fit.fitreshigh)
			p = Popen(['e/fitin', 'G0'],bufsize=0,stdin=PIPE,stdout=PIPE,stderr=PIPE,universal_newlines=True)
			print 'Automatic Sampling Started'
			profile = open('%s.%s.profile'%(axe,mol.un),'w')
			p.stdin.write('#%2d \n'%mol.id)
			l = p.stdout.read(31)
			print l
			if axe == 'X':
				v=(1,0,0)
			elif axe == 'Y':
				v=(0,1,0)
			elif axe == 'Z':
				v=(0,0,1)
			else :
				print 'strange axe ...'
				return
			mat=mol.acteur.GetMatrix()
			mat3x3 = sym.ExtractRotMat(mat)
			R = array(mat3x3)
			ivect = dot(R,v)
			(sa,sb,sg) = sym.R2Eul(mat3x3)
			angle = anglemin
			while angle < anglemax:
				gfx.sym.render=0
				mol.acteur.SetOrientation(0,0,0)
				RotaEuler(mol.acteur,sa,sb,sg)
				gfx.sym.render=1
				mol.acteur.RotateWXYZ(angle,ivect[0],ivect[1],ivect[2])
				cc = fitin_io(p,mol)
				try:
					print '%f  %f'%(angle,cc)
				except TypeError :
					MB.showwarning('Error','cc_profile : Interactive Fitting encountered a problem; go to help')
					return
				profile.write('%f  %f \n'%(angle,cc))
				angle+=step
			profile.close()
			gfx.sym.render=0
			mol.acteur.SetOrientation(0,0,0)
			RotaEuler(mol.acteur,sa,sb,sg)
			gfx.renwin.Render()
			gfx.sym.render=1
			try :
				p.stdin.write('\n')
			except :
				MB.showwarning('Error','cc_profile : Interactive Fitting terminated badly')
			print 'Automatic Sampling Stopped'
			system("gnuplot -persist <<ENDOF\n plot '%s' with line\nENDOF"%'%s.%s.profile'%(axe,mol.un))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def fitin_io(p,mol):
	mat=mol.acteur.GetMatrix()
	(a,b,g) = sym.R2Eul(sym.ExtractRotMat(mat))
	(x,y,z) = [mat.GetElement(0,3),mat.GetElement(1,3),mat.GetElement(2,3)]
	frame = '%.3f %.3f %.3f %.3f %.3f %.3f \n'%(a,b,g,x,y,z)
	try :
		p.stdin.write(frame)
		l = p.stdout.read(77)
		return float(l.split()[-1])
	except :
		MB.showwarning('Error','Problem in function ccprofile, process fitin must be dead ... RIP')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              MISC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def input_pv_file(gfx,f):
		chdir(gfx.workdir)
		try :
			pvf = open(f,'r')
		except IOError:
			MB.showwarning('Info','Problem during fitting')
			return
		molid=gfx.mol[0].id
		pvs=[]
		for ligne in pvf:
			if ligne =='':
				MB.showwarning('Info','Problem during fitting')
				return
			l=ligne.split()
			pvs+=[[molid,float(l[2]),float(l[3]),float(l[4]),float(l[5]),float(l[6]),float(l[7])]]
			molid=next_mol_id(gfx,molid) #return none if error
		for mol in gfx.mol:
			for pv in pvs:
				if mol.id==pv[0]:
					if gfx.sym != None:
						gfx.sym.render=0
					mol.acteur.SetOrientation(0,0,0)
					RotaEuler(mol.acteur,pv[1],pv[2],pv[3]) #on defini la partie rotationelle de la transformation
					mol.acteur.SetPosition(pv[4],pv[5],pv[6])
		if gfx.sym != None:
			gfx.sym.render=1
		gfx.renwin.Render()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def next_mol_id(gfx,ident):
	for i in range(len(gfx.mol)-1):
		if gfx.mol[i].id==ident:
			 return	gfx.mol[i+1].id
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def RotaEuler(acteur,alpha,beta,gamma):
		acteur.RotateWXYZ(gamma,0,0,1)
		acteur.RotateWXYZ(beta,0,1,0)
		acteur.RotateWXYZ(alpha,0,0,1)
