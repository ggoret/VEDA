# -*- coding: utf-8 -*-
# Interface Module of VEDA
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
try :
	import ttk
except :
	try :
		import tkinter.ttk as ttk
	except :
		print 'VEDA : Error : ttk is not installed or not included into Tkinter.'
import vtk
from os import chdir
import Map,gfx,mod,Nma,uro,sym
import tkFileDialog
import tkColorChooser
from ScrolledText import ScrolledText
import tkMessageBox as MB
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def quit(frame,gfx):
	if gfx.ifit != None:
		gfx.ifit.terminate_process()
	frame.quit()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Definitions de frames
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def create_frame(root):
	frame = tk.Frame(root,width=107)
	frame.pack( fill=tk.Y ,side=tk.RIGHT)
	return frame
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def create_help_frame(gfx):
	helpframe = tk.Toplevel(width=800,height=600)
	helpframe.title('Help')
        lf = tk.LabelFrame(helpframe, text= 'Help about interactor keys')
        lf.pack (side='top', fill='both', expand=1)
        lt = ScrolledText(lf,height=20,background='beige')
        lt.pack (side='top', fill='both', expand=1)
        lt.delete("0.0","end - 1 char")
        lt.insert('end',open(gfx.vedadir + '/doc/intkey.txt','r').read())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Recherche de fichier
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def browse():
	fname = tkFileDialog.askopenfilename(title="Open * file",filetypes=[("All", "*.*")])
	return fname
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def browsemol():
	fname = tkFileDialog.askopenfilename(title="Open pdb file",filetypes=[("PDB files", "*.pdb"),("Map Files", "*.ezd"),("All", "*.*")])
	return fname
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def browsemap():
	fname = tkFileDialog.askopenfilename(title="Open map file",filetypes=[("Map File", "*.ezd"),("Map File", "*.vtk"),("All", "*.*")])
	return fname
#===============================================================================
#                              CC Bar
#===============================================================================
class CCbar(vtk.vtkAxisActor2D):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self):
        self.GetProperty().SetLineWidth(15)
        self.TickVisibilityOff()
        self.LabelVisibilityOff()
        self.x=0.10
        self.y=0.15
        self.SetPosition(self.x,self.y)
        self.SetPosition2(self.x,self.y)
        self.SetTitle("")
        self.sb=CCref(self.x,self.y)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def update(self,cc):
        self.SetTitle(str(cc))
        if cc!='CC':
		try:
            		self.SetPosition2(self.x,cc/200.+self.y)
        	except:
            		print "Error in cc"
	else :
		self.SetPosition2(self.x,0+self.y)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def render(self,ren):
        ren.AddActor2D(self)
        for i in range(5):
            ren.AddActor2D(self.sb.part[i])
        self.hide()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def show(self):
        self.VisibilityOn()
        for i in range(5):
            self.sb.part[i].VisibilityOn()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def hide(self):
        self.VisibilityOff()
        for i in range(5):
            self.sb.part[i].VisibilityOff()
#===============================================================================
class CCref:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,x,y):
        self.part=[]
        for i in range(5):
            self.part.append(vtk.vtkAxisActor2D())
            self.part[i].GetProperty().SetColor(0.5,0.5,0.5)
            self.part[i].GetProperty().SetLineWidth(4)
            self.part[i].SetPosition(x,y+i*0.1)
            self.part[i].SetPosition2(x,y+.1+i*0.1)
            self.part[i].TickVisibilityOff()
            self.part[i].LabelVisibilityOff()
#===============================================================================
#                              Status Bar
#===============================================================================
class StatusBar():
  	def __init__(self,root):
		self.label=tk.Label(root,bd=1,relief=tk.SUNKEN,anchor=tk.W)
		self.label.pack(side='bottom',fill=tk.X)
  	def set(self,format,*args):
		self.label.config(background='red',text=format % args)
		self.label.update_idletasks()
  	def clear(self):
		self.label.config(text="",background='white')
		self.label.update_idletasks()
#===============================================================================
#                              Menu
#===============================================================================
class itf():
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def __init__(self,root,gfx):
		self.root=root
		self.mollist=None
		self.mollist2=None
		self.mollist3=None
		self.mollist4=None
		self.modlist=None
		self.modlist2=None
		self.modlist3=None
		self.nmmodlist=None
		self.nmmollist=None
		self.status=StatusBar(root)
		self.ccbar=CCbar()
		self.ccbar.render(gfx.renderer)
		self.debug=0
		self.nma=0
		self.embshf=None
		self.mapwizopen=0
		self.shellwizopen=0
		self.modwizopen=0
		self.molwizopen=0
		self.viscropopen=0
		self.setupwizopen=0
		self.fitwizopen=0
		self.ifitwizopen=0
		self.displaymapopen=0
		self.displaymolopen=0
		self.displayconstopen=0
		self.ccprofwizopen=0
		self.bookkeepingopen=0
		self.ncsrmswizopen=0
		self.magnumwizopen=0
		self.nmwizopen=0
		self.nmvieweropen=0
		self.nmacceptopen=0
		self.replacebut=None
		self.acceptbut=None
		self.varloop = tk.IntVar()
		self.varconformer = tk.IntVar()
		self.confscroll=None
		self.varconformerSitu=tk.IntVar()
		self.confscrollSitu = None
		self.mapscaleentry=None #entry
		self.dzentry=None #entry
		self.svg=None #StringVar
		self.fitscaleentry=None #entry
		self.centry=None #entry
		self.slabvar=tk.IntVar()
		self.zvar=tk.IntVar()
		self.varslab = None
		self.varfog=None
		self.fogdivar=tk.IntVar()
		self.fogdevar=tk.IntVar()
		self.butslab= None
		self.txtslab = None
		self.txtz = None
		self.slab = None
		self.z = None
		self.fitreshighentry=None
		self.fitreslowentry=None
		self.imodebut=None
		gfx.itf = self
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def init_menu(self,root,gfx):
		root.title("VEDA : Visual Environment for Docking Algorithms")
		#creation de l'icone
		root.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
		#Creation de la barre de menu:
		menu = tk.Menu(root)
		data = tk.Menu(menu , tearoff=0) #M#
		menu.add_cascade(label="File",menu=data)
		loadmap = tk.Menu(menu , tearoff=0) #m#
		loadmod = tk.Menu(menu , tearoff=0) #m#
		#resolution = tk.Menu(menu , tearoff=0) #m#
		data.add_command(label="Target Map", command=lambda:self.load_map_wiz(gfx))
		data.add_command(label="Models", command=lambda:self.load_models_wiz(gfx))
		data.add('separator')
		save = tk.Menu(menu, tearoff=0)
		data.add_cascade(label="Save",menu=save,state = 'active')
		save.add_command(label="Molecule", command=lambda :mod.save_mol_as(root,self.mollist,gfx))
		save.add_command(label="Constellation", command=lambda :try_save_constellation_as(self.mollist))
		save.add_command(label="Snapshot", command=lambda :gfx.snapshot())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		def try_save_constellation_as(mollist):
			if gfx.sym!=None:
				gfx.sym.save_constellation_as(mollist)
			else :
				MB.showwarning('Info','No constellation loaded')
				return
		data.add('separator')
		data.add_command(label="Exit", command=lambda:quit(root,gfx))
		molecules = tk.Menu(menu , tearoff=0) #M#
		menu.add_cascade(label="Assemblage",menu=molecules)
		molecules.add_command(label="Molecules", command=lambda:self.Manage_molecules_wiz(gfx))
		molecules.add_command(label="BookKeeping", command=lambda:self.bookkeeping(gfx))
		fitting = tk.Menu(menu , tearoff=0) #M#
		fitops = tk.Menu(menu , tearoff=0) #m#
		menu.add_cascade(label="Fitting",menu=fitting)
		fitting.add_command(label="Setup", command=lambda:self.setup_wiz(gfx))
		fitting.add('separator')
		fitting.add_command(label="Interactive Fitting", command=lambda : self.int_fit_wiz(gfx))
		fitting.add_command(label="URO Refinement", command=lambda : self.fit_wiz(gfx))
		
		fitting
		#if self.nma:
		tools = tk.Menu(menu, tearoff=0) #M#
		menu.add_cascade(label="Tools",menu=tools,state = 'active')
		tools.add_command(label="Magnification Finder", command=lambda : self.magnum_wiz(gfx))
		tools.add_command(label="RMS Threshold", command=lambda : self.ncs_rms_wiz(gfx))
		tools.add_command(label="CC profiler", command=lambda : self.cc_prof_wiz(gfx))
		tools.add('separator')
		tools.add_command(label="Normal Modes Analyses", command=lambda : self.normal_modes_wiz(gfx))
		"""
		nma.add_command(label="Build Hessian Matrix", command=lambda :Nma.pdbmat(gfx,self.mollist))
		nma.add_command(label="Diagonalize Hessian", command=lambda :Nma.diaghessian(gfx))
		nma.add_command(label="Build Nma Shifting", command=lambda :create_nma_frame(gfx))
		nma.add_command(label="Display Structure", command=lambda :Nma.display_nma(gfx))
		"""
		display = tk.Menu(menu, tearoff=0) #M#
		menu.add_cascade(label="Display",menu=display,state = 'active')
		rendumap = tk.Menu(menu, tearoff=0) #m#
		rendumol = tk.Menu(menu, tearoff=0) #m#
		#renducam = tk.Menu(menu, tearoff=0) #m#
		rendumisc = tk.Menu(menu, tearoff=0) #m#
		display.add_cascade(label="Map",menu=rendumap,state='active')
		display.add_cascade(label="Molecules",menu=rendumol,state='active')
		#display.add_cascade(label="Camera",menu=renducam,state='active')
		display.add_cascade(label="Miscelaneous",menu=rendumisc,state='active')
		rendumisc.add_command(label="Axes On/Off", command=lambda :gfx.setaxes())
		rendumisc.add_command(label="Fog On/Off", command=lambda :gfx.disable_fog())
		rendumisc.add_command(label="Visual Crop", command=lambda:self.crop_wiz(gfx))
		rendumol.add_command(label="Properties", command=lambda :self.display_mol_wiz(gfx))
		rendumap.add_command(label="Map On/Off", command=lambda :Map.hidemap(gfx))
		rendumap.add_command(label="Box On/Off", command=lambda :Map.hidemapbox(gfx))
		rendumap.add_command(label="Sym Support On/Off", command=lambda : sym.hide_show_symsupport(gfx))
		rendumap.add_command(label="Properties", command=lambda :self.display_map_wiz(gfx))
		'''
		renducam.add_command(label="Parallel Projection On/Off", command=lambda :gfx.parallel_proj_onoff())
		renducam.add_command(label="Auto Adjust Clipping On/Off", command=lambda :gfx.auto_adjust_clipping_onoff())
		renducam.add_command(label="Slab Mode On/Off", command=lambda :gfx.slab_mode_onoff())
		'''
		help = tk.Menu(menu, tearoff=0) #M#
		menu.add_cascade(label="Help",menu=help)
		help.add_command(label="Interactor Keys", command=lambda :create_help_frame(gfx))
		help.add_command(label="Display Version", command=lambda :opening())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		def opening():
			gfx.set_open_text_on_off()
			gfx.set_open_img_on_off()
		dbg = tk.Menu(menu , tearoff=0) #M#
		menu.add_cascade(label="Debug",menu=dbg)
		dbg.add_command(label="Embedded Python shell", command=lambda:gfx.embedded_shell_wiz())
		#dbg.add_command(label="Load Pdb Nma models", command=lambda : Nma.load_models_nma(gfx,menu,"") ) #""= no file = browse
		root.config(menu=menu)
		root.protocol("WM_DELETE_WINDOW", lambda:root.quit())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              listebox
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def init_mod_list(self,frame,gfx):
		txtlb = tk.Label(frame, text='Mol List :')
		txtlb.place(x=5,y=2)
		scrollbar = tk.Scrollbar(frame, orient=tk.VERTICAL)
		self.mollist = tk.Listbox(frame,selectmode=tk.SINGLE ,width=10, height=15,background='black',foreground='white')
		scrollbar.config(command=self.mollist.yview)
		self.mollist.config(yscrollcommand=scrollbar.set)
		scrollbar.place(x= 90, y=20,width=14,height=243)
		self.mollist.place(x=5,y=20)
		butlctmol = tk.Button(frame,text="Locate Mol",width=9,height=1,command=lambda :mod.locate_independent_mol(gfx))
		butlctmol.place(x= 5, y=270)
		butinvcts = tk.Button(frame,text="Const On/Off",width=9,height=1,command=lambda :mod.invert_constellation(gfx))
		butinvcts.place(x= 5, y=300)
		#butrescam = tk.Button(frame,text="Reset Cam",width=9,height=1,command=lambda :gfx.reset_view())
		#butrescam.place(x= 5, y=330)
		cam=tk.PhotoImage(file=gfx.vedadir + 'doc/cam.gif')
		self.imodebut = tk.Button(frame, width=50, height=50, image=cam,command=gfx.switch)
		self.imodebut.image=cam
		self.imodebut.pack(side=tk.BOTTOM)
		imodelb = tk.Label(frame, text=' Interactor Mode ')
		imodelb.pack(side=tk.BOTTOM)
		self.rightframe = frame
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def init_slab_stuff(self,gfx,s,z,sm,zm):
		if self.varslab == None:
			self.varslab = tk.IntVar(self.rightframe)
			self.varslab.set(0)
		if self.butslab != None :
			self.butslab.place_forget()
			self.txtslab.place_forget()
			self.txtz.place_forget()
			self.slab.place_forget()
			self.z.place_forget()
		self.butslab= tk.Checkbutton(self.rightframe, text="Slab Mode",variable=self.varslab,onvalue=1,command=lambda :gfx.slab_mode_onoff(self.varslab))
		self.butslab.place(x= 5, y=330)
		self.slabvar.set(s)
		self.zvar.set(z)
		self.txtslab = tk.Label(self.rightframe, text='Slab')
		self.txtslab.place(x=20,y=360)
		self.txtz = tk.Label(self.rightframe, text='Z Pos')
		self.txtz.place(x=65,y=360)
		self.slab = tk.Scale(self.rightframe,command=lambda e,g=gfx : g.slabchange(e),variable=self.slabvar,from_=1,to=sm,length=200,showvalue='yes',orient='vertical')
		self.slab.place(x=0,y=380)
		self.z = tk.Scale(self.rightframe,command=lambda e,g=gfx : g.zchange(e),variable=self.zvar,from_=-zm, to=zm,length=200,showvalue='yes',orient='vertical')
		self.z.place(x=50,y=380)
		if self.varslab.get()==0:
			self.slab.configure(state=tk.DISABLED)
			self.z.configure(state=tk.DISABLED)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def init_fog_stuff(self,gfx,s,z,sm,zm):
		self.varfog = tk.IntVar(self.rightframe)
		self.varfog.set(gfx.fogmode)
		self.butfog= tk.Checkbutton(self.rightframe, text="Fog",variable=self.varfog,onvalue=1,command=lambda :gfx.fog_on_off())
		self.butfog.place(x= 5, y=630)
		self.txtfogdi = tk.Label(self.rightframe, text='Dist.')
		self.txtfogdi.place(x=20,y=650)
		self.txtfogde = tk.Label(self.rightframe, text='Density')
		self.txtfogde.place(x=50,y=650)
		self.fogdivar.set(0)
		self.dist = tk.Scale(self.rightframe,command=lambda e,g=gfx : g.update_fog_dist(e),variable=self.fogdivar,from_=0,to=zm,length=125,showvalue='yes',orient='vertical')
		self.dist.place(x=0,y=680)
		self.density = tk.Scale(self.rightframe,command=lambda e,g=gfx : g.update_fog_density(e),variable=self.fogdevar,from_=100, to=0,length=125,showvalue='yes',orient='vertical')
		self.density.place(x=50,y=680)
		if self.varfog.get()==0:
			self.dist.configure(state=tk.DISABLED)
			self.density.configure(state=tk.DISABLED)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Message Box
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def showinfo(self,title,message):
		message = '\n     ' + message + '      \n'
		frame = tk.Toplevel(width=350,height=150)
		frame.title(title)
		frame.transient(self.root)
		label = tk.Label(frame,text=message,font=("bold"))
		label.pack()
		but = tk.Button(frame,width=10,text='Ok', command=lambda :frame.destroy())
		but.pack()
		frame.protocol("WM_DELETE_WINDOW", lambda:frame.destroy())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Display Map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def display_map_wiz(self,gfx):
		if not self.displaymapopen:
			self.displaymapopen=1
			frame = tk.Toplevel(width=250,height=145)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.title('Display Map')
			frame.transient(self.root)
			#initial value
			butcolor = tk.Button(frame,bg=mod.vtk2tkhex_color((0,0.5,0.75)), text='Select Color',width=10,height=1, command=lambda:Map.color_map(gfx,butcolor))
			butcolor.place( x=130,y=20)
			if gfx.map!=[]:
				butcolor.config(bg=mod.vtk2tkhex_color(gfx.map[0].color))
			txtne = tk.Label(frame, text='Rendering Type :')
			txtne.place( x=130,y=55)
			varrt = tk.StringVar(frame)
			varrt.set("Wireframe") #initial value
			rendom=tk.OptionMenu(frame,varrt,'Surface', 'Wireframe','Points')
			rendom.place(x=130,y=74)
			rendom.config(bg='beige')
			if gfx.map!=[]:
				varrt.set(gfx.map[0].rendtype)
			txtcl = tk.Label(frame, text='Contour Level :')
			txtcl.place( x=4,y=4)
			clentry = tk.Entry(frame, width=12,background='beige')
			clentry.place(  x=4,y=24)
			if gfx.map!=[]:
				clentry.delete(0,tk.END)
				clentry.insert(tk.END,str(gfx.map[0].isov))
			else :
				clentry.insert(tk.END,"1.0")
			txtcl = tk.Label(frame, text='Opacity :')
			txtcl.place( x=4,y=55)
			opentry = tk.Entry(frame, width=12,background='beige')
			opentry.place(  x=4,y=75)
			if gfx.map!=[]:
				opentry.delete(0,tk.END)
				opentry.insert(tk.END,str(gfx.map[0].opct))
			else :
				opentry.insert(tk.END,"1.0")
			butdisp = tk.Button(frame, text='Apply', command=lambda :Map.display_map(gfx,varrt.get(),clentry.get(),opentry.get(),mod.tkhex2vtk_color(butcolor['bg'])))
			butdisp.place( x=4,y=113)
			'''
			def bind_return(event):
				butdisp.invoke()
			frame.bind("<Return>",bind_return)
			'''
		else :
			MB.showwarning('Info','The Display map wizard is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.displaymapquit(frame,gfx))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def displaymapquit(self,frame,gfx):
		self.displaymapopen=0
		frame.destroy()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              MAP Wizard
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def load_map_wiz(self,gfx):
		if not self.mapwizopen:
			self.mapwizopen=1
			frame = tk.Toplevel(width=270,height=180)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.title('Target Map Wizard')
			frame.transient(self.root)
			frame.geometry("270x210+700+120")
			nb = ttk.Notebook(frame, width=270, height=180)
			nb.pressed_index = None
			frame1=tk.Frame(nb)
			butbrw = tk.Button(frame1, text='Browse Map-file', width=10,height=1,command=lambda :Map.assign_map(gfx,nfv,cpltnfv,self.mapscaleentry,clentry,opentry,butcolor))
			butbrw.place(  x=4,y=5)
			cpltnfv=nfv=tk.StringVar()
			nfv=tk.StringVar()
			nflabel = tk.Label(frame1,textvariable=nfv,width=12,background='beige')
			nflabel.place( x=5,y=35)
			if gfx.map!=[]:
				nfv.set(Map.extract_file_from_path(gfx.map[0].fn))
				cpltnfv.set(gfx.map[0].fn)
			else :
				nfv.set('')
				cpltnfv.set('')
			txtne = tk.Label(frame1, text='EM Magnification :')
			txtne.place( x=150,y=5)
			self.mapscaleentry = tk.Entry(frame1, width=12,background='beige')
			self.mapscaleentry.place(  x=150,y=25)
			if gfx.map!=[]:
				self.mapscaleentry.delete(0,tk.END)
				self.mapscaleentry.insert(tk.END,str(gfx.map[0].scale))
			else :
				self.mapscaleentry.insert(tk.END,"1.0")
			#initial value
			butcolor = tk.Button(frame1,bg=mod.vtk2tkhex_color((0,0.5,0.75)), text='Select Color',width=10,height=1, command=lambda:Map.color_map(gfx,butcolor))
			butcolor.place( x=4,y=60)
			if gfx.map!=[]:
				butcolor.config(bg=mod.vtk2tkhex_color(gfx.map[0].color))
			txtne = tk.Label(frame1, text='Rendering Type :')
			txtne.place( x=5,y=90)
			varrt = tk.StringVar(frame1)
			varrt.set("Wireframe") #initial value
			rendom=tk.OptionMenu(frame1,varrt,'Surface', 'Wireframe','Points')
			rendom.place(x=4,y=110)
			rendom.config(bg='beige')
			if gfx.map!=[]:
				varrt.set(gfx.map[0].rendtype)
			txtcl = tk.Label(frame1, text='Contour Level :')
			txtcl.place( x=150,y=55)
			clentry = tk.Entry(frame1, width=12,background='beige')
			clentry.place(  x=150,y=75)
			if gfx.map!=[]:
				clentry.delete(0,tk.END)
				clentry.insert(tk.END,str(gfx.map[0].isov))
			else :
				clentry.insert(tk.END,"1.0")
			txtcl = tk.Label(frame1, text='Opacity :')
			txtcl.place( x=150,y=106)
			opentry = tk.Entry(frame1, width=12,background='beige')
			opentry.place(  x=150,y=126)
			if gfx.map!=[]:
				opentry.delete(0,tk.END)
				opentry.insert(tk.END,str(gfx.map[0].opct))
			else :
				opentry.insert(tk.END,"1.0")
			smooth = tk.StringVar(frame1)
			smooth.set(1)
			smoothing= tk.Checkbutton(frame1, text="Smoothing",variable=smooth,onvalue=1)
			smoothing.place(x=70,y=150)
			if gfx.map!=[]:
				smooth.set(gfx.map[0].issmooth)
			vardeci = tk.StringVar(frame1)
			vardeci.set(0)
			cbig= tk.Checkbutton(frame1, text="Decimation",variable=vardeci,onvalue=1)
			cbig.place(x=160,y=150)
			if gfx.map!=[]:
				vardeci.set(gfx.map[0].isdeci)
			butdisp = tk.Button(frame1, text='Load', command=lambda :Map.load_map(gfx,cpltnfv.get(),self.root,self.status,float(self.mapscaleentry.get()),varrt.get(),float(clentry.get()),float(opentry.get()),cropentry,nfv,vardeci.get(),smooth.get(),mod.tkhex2vtk_color(butcolor['bg']),caller='map'))
			butdisp.place( x=4,y=145)
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	#                        SYMMETRY                                  #
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
			frame2 = tk.Frame(nb)
			txtsf = tk.Label(frame2, text='Select Symmetry :')
			txtsf.place( x=4,y=4)
			varsym = tk.StringVar(frame2)
			if gfx.ps!=None:
				varsym.set(gfx.ps.solidtype)
			else :
				varsym.set("P1") #initial value
			symom=tk.OptionMenu(frame2,varsym,'Icosahedral','Octahedral','Tetrahedral','Helicoidal', 'Cn', 'Dn', 'P1','User Defined')
			symom.place(x=4,y=22)
			#automatic sym settings
			symom['menu'].entryconfigure('Icosahedral',command=lambda:display_sym_option(self.sgv,2,'Icosahedral'))
			symom['menu'].entryconfigure('Octahedral',command=lambda:display_sym_option(self.sgv,2,'Octahedral'))
			symom['menu'].entryconfigure('Tetrahedral',command=lambda:display_sym_option(self.sgv,2,'Tetrahedral'))
			symom['menu'].entryconfigure('Helicoidal',command=lambda:display_sym_option(self.sgv,2,'Helicoidal'))
			symom['menu'].entryconfigure('Dn',command=lambda:display_sym_option(self.sgv,2,'Dn'))
			symom['menu'].entryconfigure('Cn',command=lambda:display_sym_option(self.sgv,2,'Cn'))
			symom['menu'].entryconfigure('P1',command=lambda:display_sym_option(self.sgv,2,'P1'))
			symom['menu'].entryconfigure('User Defined',command=lambda:display_sym_option(self.sgv,2,'User Defined'))
			suprad = tk.Label(frame2, text='Support Radius :')
			self.sgv=tk.StringVar() #hint
			srentry = tk.Entry(frame2,textvariable=self.sgv, width=12,background='beige')
			if gfx.ps!=None:
				if gfx.ps.radius!=None:
					srentry.delete(0,tk.END)
					srentry.insert(tk.END,str(gfx.ps.radius))
				else:
					srentry.insert(tk.END,"200")
			else :
				srentry.insert(tk.END,"200")
			butaply = tk.Button(frame2, text='Load', command=lambda : load_sym()) #sr = support radius
			butaply.place(x=4,y=145)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def load_sym():
				sym.ps(gfx,varsym.get(),srentry,varori.get(),cTUentry,deltaentry,varttub.get(),dnnentry.get(),cnnentry.get(),varaxe.get(),symfilevar,phizeroentry.get(),caller='sym')
				if varsym.get()=='Helicoidal':
					varttub.set('Elm-Hel')
					gfx.ps.enttype='Elm-Hel'
					display_sym_option(self.sgv,2,'Helicoidal')
			#~~~~~~set_ori~~~~~~~~~~~~~~~~~~~~~~~~~~#
			varori = tk.StringVar(frame2) #ori for platonic
			if gfx.ps!=None :
				varori.set(gfx.ps.ori)
			else :
				varori.set("") #initial value
			#~~~~~~ico~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
			symopico = tk.Label(frame2, text='Ico Settings :')
			symopico.place( x=4,y=55)
			symopico.place_forget()
			icoorilst=['5Z2Y.1','5Z2Y.2','5Z2X.1','5Z2X.2','5Y2Z.1','5Y2Z.2','5Y2X.1','5Y2X.2','5X2Z.1','5X2Z.2','5X2Y.1','5X2Y.2','3Z2Y.1','3Z2Y.2','3Z2X.1','3Z2X.2','3Y2Z.1','3Y2Z.2','3Y2X.1','3Y2X.2','3X2Z.1','3X2Z.2','3X2Y.1','3X2Y.2','2Z2Y.1','2Z2Y.2']
			icoom=	tk.OptionMenu(frame2,varori,*icoorilst)
			icoom.place(x=4,y=75)
			icoom.place_forget()
			#~~~~~~octa~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
			symopocta = tk.Label(frame2, text='Octa Settings :')
			symopocta.place( x=4,y=55)
			symopocta.place_forget()
			octaorilst=['4Z4Y','4Z2Y','4Y2Z','4X2Z','3Z2Y.1','3Z2Y.2','3Z2X.1','3Z2X.2','3Y2Z.1','3Y2Z.2','3Y2X.1','3Y2X.2','3X2Z.1','3X2Z.2','3X2Y.1','3X2Y.2']
			octaom=	tk.OptionMenu(frame2,varori,*octaorilst)
			octaom.place(x=4,y=75)
			octaom.place_forget()
			#~~~~~~tetra~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
			symoptetra = tk.Label(frame2, text='Tetra Settings :')
			symoptetra.place( x=4,y=55)
			symoptetra.place_forget()
			tetraorilst = ['2Z2Y2X']
			tetraom=tk.OptionMenu(frame2,varori,*tetraorilst)
			tetraom.place(x=4,y=75)
			tetraom.place_forget()
			#~~~~~~tubes~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
			varttub = tk.StringVar(frame2)
			if gfx.ps !=None:
				if gfx.ps.enttype!=None:
					varttub.set(gfx.ps.enttype)
				else :
					varttub.set('')
			else :
				varttub.set('') #initial value
			symoptub = tk.Label(frame2, text='Helical Settings :')
			ttubom=tk.OptionMenu(frame2,varttub,'Elm-Hel','cTU')
			ttubom.place(x=4,y=75)
			ttubom['menu'].entryconfigure('Elm-Hel',command=lambda:display_sym_tube_option('Elm-Hel'))
			ttubom['menu'].entryconfigure('cTU',command=lambda:display_sym_tube_option('cTU'))
			ttubom.place_forget()
			phizerolab = tk.Label(frame2, text='Φ Start :')
			phizeroentry = tk.Entry(frame2,width=6, background='beige')
			if gfx.ps!=None and gfx.ps.solidtype == 'Helicoidal':
				phizeroentry.delete(0,tk.END)
				phizeroentry.insert(tk.END,str(gfx.ps.phizero))
			else :
				phizeroentry.insert(tk.END,str(0))
			clab = tk.Label(frame2, text='c :')
			self.centry = tk.Entry(frame2,width=7, background='beige')
			Tlab = tk.Label(frame2, text='T :')
			Tentry = tk.Entry(frame2,width=5, background='beige')
			Ulab = tk.Label(frame2, text='U :')
			Uentry = tk.Entry(frame2,width=5, background='beige')
			dphilab = tk.Label(frame2, text='ΔΦ :')
			dphientry = tk.Entry(frame2,width=7, background='beige')
			dzlab = tk.Label(frame2, text='Δz :')
			self.dzentry = tk.Entry(frame2,width=7, background='beige')
			slab = tk.Label(frame2, text='#s :')
			sentry = tk.Entry(frame2,width=3, background='beige')
			if gfx.ps!=None and gfx.ps.solidtype == 'Helicoidal':
				sentry.delete(0,tk.END)
				sentry.insert(tk.END,str(gfx.ps.s))
				self.dzentry.delete(0,tk.END)
				self.dzentry.insert(tk.END,'%.3f'%gfx.ps.dz)
				dphientry.delete(0,tk.END)
				dphientry.insert(tk.END,'%.3f'%gfx.ps.dphi)
				Uentry.delete(0,tk.END)
				Uentry.insert(tk.END,str(gfx.ps.U))
				Tentry.delete(0,tk.END)
				Tentry.insert(tk.END,str(gfx.ps.T))
				self.centry.delete(0,tk.END)
				self.centry.insert(tk.END,str(gfx.ps.c))
			cTUentry = (self.centry,Tentry,Uentry)
			deltaentry=(dphientry,self.dzentry,sentry)
			#~~~~~~Dn~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
			symlabdn = tk.Label(frame2, text='Dn Settings :')
			dnnlab = tk.Label(frame2, text='n :')
			dnnentry = tk.Entry(frame2,width=5, background='beige')
			if gfx.ps!=None and gfx.ps.solidtype == 'Dn':
				dnnentry.delete(0,tk.END)
				dnnentry.insert(tk.END,str(gfx.ps.dnn))
			#~~~~~~Cn~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
			symlabcn = tk.Label(frame2, text='Cn Settings :')
			cnnlab = tk.Label(frame2, text='n :')
			cnnentry = tk.Entry(frame2,width=5, background='beige')
			phizerolab = tk.Label(frame2, text='Φ Start :')
			phizeroentry = tk.Entry(frame2,width=6, background='beige')
			varaxe = tk.StringVar(frame2)
			axelab = tk.Label(frame2, text='Axe :')
			axeom=	tk.OptionMenu(frame2,varaxe,'X','Y', 'Z')
			if gfx.ps!=None and gfx.ps.solidtype == 'Cn':
				phizeroentry.delete(0,tk.END)
				phizeroentry.insert(tk.END,str(gfx.ps.phizero))
				cnnentry.delete(0,tk.END)
				cnnentry.insert(tk.END,str(gfx.ps.cnn))
				varaxe.set(gfx.ps.axe)
			else :
				phizeroentry.insert(tk.END,str(0))
				varaxe.set('Z') #initial value
			#~~~~~~User Defined~~~~~~~~~~~~~~~~~~~~~#
			butbrwsym = tk.Button(frame2, text='Browse Sym-file', width=10,height=1,command=lambda :sym.browsesymfile(gfx,symfilevar))
			symfilevar=tk.StringVar()
			symfilelab = tk.Label(frame2,textvariable=symfilevar,anchor=tk.E,width=16,background='beige')
			if gfx.ps!=None and gfx.ps.solidtype == 'User Defined':
				symfilevar.set(gfx.ps.symfile)
			else :
				symfilevar.set("")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              dynamic display
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def display_sym_option(sgv,mode,stype):
				if mode == 2:
					varsym.set(stype)
				if varsym.get()=='P1':
					forget_all()
				if varsym.get()=='Icosahedral':
					#hint
					if gfx.map != []:
						if gfx.map[0].box!= None:
							sgv.set('%.2f'%abs(max(gfx.map[0].box.GetBounds())))
					if gfx.ps != None:
						if gfx.ps.radius != None and gfx.map!=[]:
							sgv.set('%.2f'%(gfx.ps.radius *  gfx.map[0].scale))
					#forget
					forget_all()
					#remember
					symopico.place( x=4,y=55)
					icoom.place(x=4,y=75)
					suprad.place( x=150,y=5)
					srentry.place(x=150,y=25)
					if varori.get() not in icoorilst:
						varori.set('')
				if varsym.get()=='Octahedral':
					#hint
					if gfx.map != []:
						if gfx.map[0].box!= None:
							sgv.set('%.2f'%abs(max(gfx.map[0].box.GetBounds())))
					if gfx.ps != None:
						if gfx.ps.radius != None and gfx.map!=[]:
							sgv.set('%.2f'%(gfx.ps.radius *  gfx.map[0].scale))
					#forget
					forget_all()
					#remember
					symopocta.place( x=4,y=55)
					octaom.place(x=4,y=75)
					suprad.place( x=150,y=5)
					srentry.place(x=150,y=25)
					if varori.get() not in octaorilst:
						varori.set('')
				if varsym.get()=='Tetrahedral':
					#hint
					if gfx.map != []:
						if gfx.map[0].box!= None:
							sgv.set('%.2f'%abs(max(gfx.map[0].box.GetBounds())))
					if gfx.ps != None:
						if gfx.ps.radius != None and gfx.map!=[]:
							sgv.set('%.2f'%(gfx.ps.radius *  gfx.map[0].scale))
					#forget
					forget_all()
					#remember
					symoptetra.place( x=4,y=55)
					tetraom.place(x=4,y=75)
					suprad.place( x=150,y=5)
					srentry.place(x=150,y=25)
					if varori.get() not in tetraorilst:
						varori.set('')
				if varsym.get()=='Helicoidal':
					#hint
					if gfx.map != []:
						if gfx.map[0].box!= None:
							sgv.set('%.2f'%(abs(max(gfx.map[0].box.GetBounds())/2)))
					if gfx.ps != None:
						if gfx.ps.radius != None and gfx.map!=[]:
							sgv.set('%.2f'%(gfx.ps.radius *  gfx.map[0].scale))
					#forget
					forget_all()
					#remember
					symoptub.place( x=4,y=55)
					ttubom.place(x=4,y=75)
					display_sym_tube_option(varttub.get())
					suprad.place( x=150,y=5)
					srentry.place(x=150,y=25)
				if varsym.get()=='Dn':
					#hint
					if gfx.map != []:
						if gfx.map[0].box!= None:
							sgv.set('%.2f'%abs(max(gfx.map[0].box.GetBounds())))
					if gfx.ps != None:
						if gfx.ps.radius != None and gfx.map!=[]:
							sgv.set('%.2f'%(gfx.ps.radius *  gfx.map[0].scale))
					#forget
					forget_all()
					#remember
					symlabdn.place( x=4,y=55)
					dnnlab.place(x=4,y=75)
					dnnentry.place(x=24,y=75)
					suprad.place( x=150,y=5)
					srentry.place(x=150,y=25)
				if varsym.get()=='Cn':
					#hint
					if gfx.map != []:
						if gfx.map[0].box!= None:
							sgv.set('%.2f'%abs(max(gfx.map[0].box.GetBounds())))
					if gfx.ps != None:
						if gfx.ps.radius != None and gfx.map!=[]:
							sgv.set('%.2f'%(gfx.ps.radius *  gfx.map[0].scale))
					#delta phi
					if gfx.ps!=None :
						dphientry.delete(0,tk.END)
						if type(gfx.ps.dphi)!=str:
							dphientry.insert(tk.END,'%.3f'%gfx.ps.dphi)
						else :
							dphientry.insert(tk.END,gfx.ps.dphi)
					#forget
					forget_all()
					#remember
					phizerolab.place( x=190,y=55)
					phizeroentry.place(x=190,y=75)
					symlabcn.place( x=4,y=55)
					cnnlab.place(x=4,y=75)
					cnnentry.place(x=24,y=75)
					axelab.place(x=80,y=75)
					axeom.place(x=115,y=70)
					suprad.place( x=150,y=5)
					srentry.place(x=150,y=25)
				if varsym.get()=='User Defined':
					#forget
					forget_all()
					#remember
					butbrwsym.place( x=4,y=55)
					symfilelab.place( x=4,y=85)
				else:
					return
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def forget_all():
				symoptub.place_forget()
				phizerolab.place_forget()
				phizeroentry.place_forget()
				clab.place_forget()
				self.centry.place_forget()
				Tlab.place_forget()
				Tentry.place_forget()
				Ulab.place_forget()
				Uentry.place_forget()
				dphilab.place_forget()
				dphientry.place_forget()
				dzlab.place_forget()
				self.dzentry.place_forget()
				ttubom.place_forget()
				slab.place_forget()
				sentry.place_forget()
				symopico.place_forget()
				icoom.place_forget()
				symopocta.place_forget()
				octaom.place_forget()
				suprad.place_forget()
				srentry.place_forget()
				dnnlab.place_forget()
				dnnentry.place_forget()
				symlabdn.place_forget()
				symlabcn.place_forget()
				cnnlab.place_forget()
				cnnentry.place_forget()
				axeom.place_forget()
				axelab.place_forget()
				butbrwsym.place_forget()
				symfilelab.place_forget()
				symopocta.place_forget()
				octaom.place_forget()
				symoptetra.place_forget()
				tetraom.place_forget()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def display_sym_tube_option(enttype):
				if enttype=='cTU':
					varttub.set('cTU')
					phizerolab.place( x=190,y=55)
					phizeroentry.place(x=190,y=75)
					clab.place(x=25,y=110)
					self.centry.place(x=41,y=110)
					Tlab.place(x=115,y=110)
					Tentry.place(x=131,y=110)
					Ulab.place(x=185,y=110)
					Uentry.place(x=206,y=110)
					dphilab.place_forget()
					dphientry.place_forget()
					dzlab.place_forget()
					self.dzentry.place_forget()
					slab.place_forget()
					if gfx.ps!=None :
						sentry.delete(0,tk.END)
						sentry.insert(tk.END,'1')
						gfx.ps.s = 1
					sentry.place_forget()
				elif enttype=='Elm-Hel':
					varttub.set('Elm-Hel')
					phizerolab.place( x=190,y=55)
					phizeroentry.place(x=190,y=75)
					dphilab.place(x=5,y=110)
					dphientry.place(x=32,y=110)
					if gfx.ps!=None :
						dphientry.delete(0,tk.END)
						if type(gfx.ps.dphi)!=str:
							dphientry.insert(tk.END,'%.3f'%gfx.ps.dphi)
						else :
							dphientry.insert(tk.END,gfx.ps.dphi)
					dzlab.place(x=97,y=110)
					self.dzentry.place(x=122,y=110)
					if gfx.ps!=None:
						self.dzentry.delete(0,tk.END)
						if type(gfx.ps.dz)!=str:
							self.dzentry.insert(tk.END,'%.3f'%gfx.ps.dz)
						else :
							self.dzentry.insert(tk.END,gfx.ps.dz)
					slab.place(x=190,y=110)
					sentry.place(x=220,y=110)
					if gfx.ps!=None :
						sentry.delete(0,tk.END)
						sentry.insert(tk.END,str(gfx.ps.s))
					clab.place_forget()
					self.centry.place_forget()
					Tlab.place_forget()
					Tentry.place_forget()
					Ulab.place_forget()
					Uentry.place_forget()
			display_sym_option(self.sgv,0,'') #0 = automatique,'' = stype vide
#                              CROP
			frame3 = tk.Frame(nb)
			butshow3 = tk.Button(frame3, text='Start', command=lambda :start_crop(gfx))
			butshow3.place(  x=25,y=5)
			butstop = tk.Button(frame3, text='Stop', command=lambda : stop_crop(gfx))
			butstop.place(  x=90,y=5)
			butreset = tk.Button(frame3, text='Reset', command=lambda : reset_crop(gfx))
			butreset.place(  x=155,y=5)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def start_crop(gfx):
				try:
					butcrop.configure(state=tk.NORMAL)
					gfx.crop.show()
				except AttributeError:
					MB.showwarning('Info','Nothing to crop, load a map')
					return
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def stop_crop(gfx):
				butcrop.configure(state=tk.DISABLED)
				if gfx.crop!=None:
					gfx.crop.reset_clipping()
					gfx.crop.hide()
				else :
					if nb.index("current")==2:
						MB.showwarning('Info','Start, first')
						return
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def reset_crop(gfx):
				if gfx.crop!=None:
					gfx.crop.reset_crop_bounds()
					gfx.crop.reset_clipping()
				else :
					MB.showwarning('Info','Reset, first')
					return
			xllab = tk.Label(frame3, text='x low :')
			xllab.place( x=5,y=45)
			xlentry = tk.Entry(frame3,width=7, background='beige')
			xlentry.place(x=50,y=45)
			yllab = tk.Label(frame3, text='y low :')
			yllab.place( x=5,y=70)
			ylentry = tk.Entry(frame3,width=7, background='beige')
			ylentry.place(x=50,y=70)
			zllab = tk.Label(frame3, text='z low :')
			zllab.place( x=5,y=95)
			zlentry = tk.Entry(frame3,width=7, background='beige')
			zlentry.place(x=50,y=95)
			xulab = tk.Label(frame3, text='x up :')
			xulab.place( x=130,y=45)
			xuentry = tk.Entry(frame3,width=7, background='beige')
			xuentry.place(x=170,y=45)
			yulab = tk.Label(frame3, text='y up :')
			yulab.place( x=130,y=70)
			yuentry = tk.Entry(frame3,width=7, background='beige')
			yuentry.place(x=170,y=70)
			zulab = tk.Label(frame3, text='z up :')
			zulab.place( x=130,y=95)
			zuentry = tk.Entry(frame3,width=7, background='beige')
			zuentry.place(x=170,y=95)
			cropentry=(xlentry,xuentry,ylentry,yuentry,zlentry,zuentry)
			if gfx.crop!=None:
				xlentry.delete(0,tk.END)
				xlentry.insert(tk.END,str(gfx.crop.entryval[0]))
				xuentry.delete(0,tk.END)
				xuentry.insert(tk.END,str(gfx.crop.entryval[1]))
				ylentry.delete(0,tk.END)
				ylentry.insert(tk.END,str(gfx.crop.entryval[2]))
				yuentry.delete(0,tk.END)
				yuentry.insert(tk.END,str(gfx.crop.entryval[3]))
				zlentry.delete(0,tk.END)
				zlentry.insert(tk.END,str(gfx.crop.entryval[4]))
				zuentry.delete(0,tk.END)
				zuentry.insert(tk.END,str(gfx.crop.entryval[5]))
				gfx.crop=Map.Crop(gfx,gfx.map[0].iso,cropentry,gfx.crop.entryval)
			butcrop = tk.Button(frame3, text='Crop Map-File', command=lambda : cropandload())
			butcrop.place(  x=50,y=140)
			butcrop.configure(state=tk.DISABLED)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def cropandload():
				if gfx.map!=[]:
					gfx.crop.crop_file()
					scale = gfx.map[0].scale
					Map.load_map(gfx,gfx.crop.crpfn,self.root,self.status,scale,varrt.get(),float(clentry.get()),float(opentry.get()),cropentry,nfv,vardeci.get(),smooth.get(),gfx.map[0].color,caller='crop')
				else:
					MB.showwarning('Info','Loading a map')
					return
#                              TABS MANAGER
			'''
			def bind_return(event):
				if nb.index("current")==0:
					butdisp.invoke()
				elif nb.index("current")==1:
					butaply.invoke()
				elif nb.index("current")==2:
					butcrop.invoke()
			frame.bind("<Return>",bind_return)
			'''
			def bind_change_to_sym():
				if nb.index("current")==0:
					stop_crop(gfx)
				if nb.index("current")==1:
					stop_crop(gfx)
					if gfx.ps != None:
						display_sym_option(self.sgv,2,gfx.ps.solidtype)
			nb.bind('<<NotebookTabChanged>>', lambda evt: bind_change_to_sym())
			nb.add(frame1, text='Map', padding=1)
			nb.add(frame2, text='Symmetry', padding=1)
			nb.add(frame3, text='Crop', padding=1)
			nb.pack(expand=1, fill='both')
		else :
			MB.showwarning('Info','The target map wizard is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.mapwizquit(frame,gfx))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def mapwizquit(self,frame,gfx):
		if gfx.crop != None:
			#gfx.crop.reset_crop_bounds()
			gfx.crop.reset_clipping()
			gfx.crop.hide()
		frame.destroy()
		self.mapwizopen=0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Model Wizard
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def load_models_wiz(self,gfx):
		if not self.modwizopen:
			self.modwizopen=1
			frame = tk.Toplevel(width=300,height=515)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.geometry("+700+120")
			frame.transient(self.root)
			frame.title('Models Wizard')
			butbrw = tk.Button(frame, text='Browse File', width=10,height=1,command=lambda :assign_and_display_mod(gfx,nfv,cpltnfv,nameentry))
			butbrw.place(  x=4,y=10)
			def display_mod_typ_om(nfv):
				modtyplab.place_forget()
				repmolom.place_forget()
				txtcl.place_forget()
				clentry.place_forget()
				if nfv.get().lower().endswith('pdb') :
					varrep.set('C alpha') #initial value
					modtyplab.place(x=4,y=113)
					repmolom.place(x=4,y=130)
					clentry.delete(0,tk.END)
					clentry.insert(tk.END,"-1.0")
				elif nfv.get().lower().endswith('ezd'): ### dans le future inclure ici les nouveaux formats pour les cartes
					varrep.set('Surface') #initial value
					txtcl.place( x=4,y=113)
					clentry.place(  x=4,y=130)
					clentry.delete(0,tk.END)
					clentry.insert(tk.END,"1.0")
			def assign_and_display_mod(gfx,nfv,cpltnfv,nameentry):
					butadd2lst.configure(state=tk.NORMAL)
					mod.assign_mod(gfx,nfv,cpltnfv,nameentry)
					display_mod_typ_om(nfv)
			cpltnfv=tk.StringVar() #contient le path complet
			nfv=tk.StringVar()
			nflabel = tk.Label(frame,textvariable=nfv,width=17,background='beige')
			nflabel.place( x=5,y=44)
			varrep = tk.StringVar(frame)
			modtyplab = tk.Label(frame, text='Select Model Type :')
			repmolom=tk.OptionMenu(frame,varrep,'All Atoms','Backbone','C alpha')
			txtcl = tk.Label(frame, text='Contour Level :')
			clentry = tk.Entry(frame, width=12,background='beige')
			txtne = tk.Label(frame, text='Enter Model Name :')
			txtne.place( x=4,y=65)
			nameentry = tk.Entry(frame, width=14,background='beige')
			nameentry.place(  x=4,y=85)
			butadd2lst = tk.Button(frame, text='Add -->',width=10,height=1, command=lambda :mod.add_mod(gfx,self.modlist,nameentry,nfv,cpltnfv,varrep,clentry.get()))
			butadd2lst.place( x=4,y=170)
			'''
			def bind_return(event):
				butadd2lst.invoke()
			frame.bind("<Return>",bind_return)
			'''
			butrm2lst = tk.Button(frame, text='<-- Remove',width=10,height=1, command=lambda :mod.remove_mod(gfx,self.modlist,nameentry,nfv))
			butrm2lst.place(x=160,y=170)
			txtlb = tk.Label(frame, text='Models List :')
			txtlb.place(x=170,y=3)
			self.modlist=tk.Listbox(frame,selectmode=tk.EXTENDED,width=12, height=8)
			scrollbar = tk.Scrollbar(frame, orient=tk.VERTICAL)
			scrollbar.config(command=self.modlist.yview)
			self.modlist.config(yscrollcommand=scrollbar.set)
			scrollbar.place(x= 272, y=26,width=14,height=130)
			self.modlist.place(x=170,y=25)
			if gfx.mod!=[]:
				for modl in gfx.mod:
					self.modlist.insert(tk.END,modl.un)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def display_info(evt):
				self.root.configure(cursor='watch')
				self.status.set('Display loading ... please wait')
				butadd2lst.configure(state=tk.DISABLED)
				modtyplab.place_forget()
				repmolom.place_forget()
				txtcl.place_forget()
				clentry.place_forget()
				for selmod in self.modlist.curselection():
					for mod in gfx.mod:
						if self.modlist.get(selmod)==mod.un:
							nfv.set(Map.extract_file_from_path(mod.rfn))
							nameentry.delete(0,tk.END)
							nameentry.insert(tk.END,mod.un)
							gfx.miniwin.display_mod(mod)
				self.status.clear()
				self.root.configure(cursor='arrow')
			self.modlist.bind("<ButtonRelease-1>",display_info)
			gfx.miniwin = mod.miniwin(frame)
			frame.geometry("300x515+700+120")
		else :
			MB.showwarning('Info','The model wizard is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.modwizquit(frame,gfx))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def modwizquit(self,frame,gfx):
		self.modwizopen=0
		del gfx.miniwin
		self.modlist=None
		frame.destroy()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              DISPLAY Molecule Wizzard
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def display_mol_wiz(self,gfx):
		if not self.displaymolopen:
			frame = tk.Toplevel(width=250,height=170)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.title('Molecules Display')
			frame.transient(self.root)
			txtsmn = tk.Label(frame, text='Select Molecule :')
			txtsmn.place( x=4,y=4)
			mollist=tk.Listbox(frame,selectmode=tk.EXTENDED,width=13, height=8) #,background='black',foreground='white')
			scrollbar = tk.Scrollbar(frame, orient=tk.VERTICAL)
			scrollbar.config(command=mollist.yview)
			mollist.config(yscrollcommand=scrollbar.set)
			scrollbar.place(x= 105, y=26,width=15,height=130)
			mollist.place(x=4,y=25)
			for mol in gfx.mol:
				mollist.insert(tk.END,mol.un)
				mollist.itemconfigure(tk.END,background=mod.vtk2tkhex_color(mol.col))
			mollist.select_set(0,tk.END)
			constlab = tk.Label(frame, text='Constellation :')
			constlab.place( x=125,y=4)
			butshowconst = tk.Button(frame, text='Show', command=lambda :mod.show_constellation(gfx,mollist))
			butshowconst.place(x=125,y=25)
			buthideconst = tk.Button(frame, text='Hide', command=lambda :mod.hide_constellation(gfx,mollist))
			buthideconst.place(x=190,y=25)
			mollab = tk.Label(frame, text='Molecules :')
			mollab.place( x=125,y=60)
			butshow = tk.Button(frame, text='Show', command=lambda :mod.show_mol(gfx,mollist))
			butshow.place( x=125,y=81)
			buthide = tk.Button(frame, text='Hide', command=lambda :mod.hide_mol(gfx,mollist))
			buthide.place( x=190,y=81)
			colorlab = tk.Label(frame, text='Color :')
			colorlab.place( x=125,y=115)
			varcolor = tk.StringVar(frame)
			butcolor = tk.OptionMenu(frame,varcolor,'red', 'green','blue','purple', 'yellow','indigo','user defined')
			butcolor.place(x=125,y=135)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def set_om_user_color():
				new_color = tkColorChooser.askcolor (title="define user color")
				if new_color[1] != None:
					    	col = mod.tk2vtk_color (new_color[0])
					    	varcolor.set('user defined')
			    			butcolor['menu'].entryconfigure('user defined',background=mod.vtk2tkhex_color(col))
						butcolor.configure(background=mod.vtk2tkhex_color(col))
						mod.display_mol(gfx,mollist,butcolor)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def set_om_color(item,col):
					varcolor.set(item)
					butcolor.configure(background=col)
					mod.display_mol(gfx,mollist,butcolor)
			butcolor['menu'].entryconfigure('red',background='#ff393c',command=lambda:set_om_color('red','#ff393c'))
			butcolor['menu'].entryconfigure('green',background='#46ff76',command=lambda:set_om_color('green','#46ff76'))
			butcolor['menu'].entryconfigure('blue',background='#0098d8',command=lambda:set_om_color('blue','#0098d8'))
			butcolor['menu'].entryconfigure('purple',background='#ff5cec',command=lambda:set_om_color('purple','#ff5cec'))
			butcolor['menu'].entryconfigure('yellow',background='#f9f341',command=lambda:set_om_color('yellow','#f9f341'))
			butcolor['menu'].entryconfigure('indigo',background='#5cf3e6',command=lambda:set_om_color('indigo','#5cf3e6'))
			butcolor['menu'].entryconfigure('user defined',command=lambda:set_om_user_color())
			self.displaymolopen=1
		else :
			MB.showwarning('Info','The Display molecule wizard is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.displaymolquit(frame,gfx))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def displaymolquit(self,frame,gfx):
		self.displaymolopen=0
		frame.destroy()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Molecule Wizzard
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def Manage_molecules_wiz(self,gfx):
		if not self.molwizopen:
			frame = tk.Toplevel(width=400,height=235)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.geometry("400x235+600+120")
			frame.title('Molecules Wizard')
			frame.transient(self.root)
			nb = ttk.Notebook(frame, width=270, height=180)
			nb.pressed_index = None
			frame1=tk.Frame(nb)
			txtsmn = tk.Label(frame1, text='Select              :')
			txtsmn.place( x=4,y=2)
			txtsmn2 = tk.Label(frame1, text='Model',font=("bold"))
			txtsmn2.place( x=45,y=0)
			self.modlist2=tk.Listbox(frame1,selectmode=tk.SINGLE ,width=12, height=5)
			scrollbar = tk.Scrollbar(frame1, orient=tk.VERTICAL)
			scrollbar.config(command=self.modlist2.yview)
			self.modlist2.config(yscrollcommand=scrollbar.set)
			scrollbar.place(x= 105, y=26,width=14,height=80)
			self.modlist2.place(x=4,y=25)
			for modo in gfx.mod:
				self.modlist2.insert(tk.END,modo.un)
			def display_mol_name(event):
				try:
					if self.modlist2.curselection()==():
						return
					for mod in gfx.mod:
						if mod.un == self.modlist2.get(tk.ACTIVE):
							nameentry.delete(0,tk.END)
							nameentry.insert(tk.END,'%s_%d'%(mod.un,mod.nmol+1))
							break
					for mod in gfx.mod:
						if self.modlist2.get(tk.ACTIVE)==mod.un:
							if mod.type=='map':
								if rendvar.get() not in ['Surface', 'Wireframe','Points']:
									rendmolom.place_forget()
									rendmapom.place_forget()
									rendvar.set("Surface")
									rendmapom.place(x=135,y=23)
								else : pass
							elif mod.type=='mol':
								if rendvar.get() not in ['Line', 'Ball','Stick','Ball & Line','Ball & Stick']:
									rendmolom.place_forget()
									rendmapom.place_forget()
									rendvar.set("Line")
									rendmolom.place(x=135,y=23)
								else : pass
							else :
								print "strange Error in gfx.mol[tk.ACTIVE].mod.type"
				except :
					pass
			def focus(event):
				self.modlist2.focus_get()
			self.modlist2.bind('<Enter>', focus)
			self.modlist2.bind_all('<ButtonRelease>', display_mol_name)
			txtne = tk.Label(frame1, text='Rendering :')
			txtne.place( x=135,y=2)
			rendvar = tk.StringVar(frame1)
			rendmolom=tk.OptionMenu(frame1,rendvar,'Line', 'Ball','Stick','Ball & Line','Ball & Stick')
			rendmapom=tk.OptionMenu(frame1,rendvar,'Surface', 'Wireframe','Points')
			colorlab = tk.Label(frame1, text='Color :')
			colorlab.place( x=135,y=56)
			varcolor = tk.StringVar(frame1)
			truecolor = tk.IntVar(frame1)
			truecolor.set(0)
			txtne2 = tk.Label(frame1, text='Enter                    Name :',)
			txtne2.place( x=4,y=118)
			txtne = tk.Label(frame1, text='Molecule',font=("bold"))
			txtne.place( x=40,y=116)
			nameentry = tk.Entry(frame1,background='beige')
			nameentry.place(x=4,y=140)
			butcolor = tk.OptionMenu(frame1,varcolor,'red', 'green','purple', 'yellow','indigo','blue','user defined','true colors')
			butcolor.place(x=135,y=76)
			def set_om_user_color():
				new_color = tkColorChooser.askcolor (title="define user color")
				if new_color[1] != None:
					    	col = mod.tk2vtk_color (new_color[0])
					    	varcolor.set('user defined')
			    			butcolor['menu'].entryconfigure('user defined',background=mod.vtk2tkhex_color(col))
						butcolor.configure(background=mod.vtk2tkhex_color(col))
						truecolor.set(0)
			def set_om_color(item,col,tc):
					varcolor.set(item)
					butcolor.configure(background=col)
					truecolor.set(tc)
			allcol =[('red','#ff393c',0),('green','#46ff76',0),('purple','#ff5cec',0),('yellow','#f9f341',0),('indigo','#5cf3e6',0),('blue','#0098d8',0)]
			butcolor['menu'].entryconfigure('red',background='#ff393c',command=lambda:set_om_color('red','#ff393c',0))
			butcolor['menu'].entryconfigure('green',background='#46ff76',command=lambda:set_om_color('green','#46ff76',0))
			butcolor['menu'].entryconfigure('blue',background='#0098d8',command=lambda:set_om_color('blue','#0098d8',0))
			butcolor['menu'].entryconfigure('purple',background='#ff5cec',command=lambda:set_om_color('purple','#ff5cec',0))
			butcolor['menu'].entryconfigure('yellow',background='#f9f341',command=lambda:set_om_color('yellow','#f9f341',0))
			butcolor['menu'].entryconfigure('indigo',background='#5cf3e6',command=lambda:set_om_color('indigo','#5cf3e6',0))
			butcolor['menu'].entryconfigure('user defined',command=lambda:set_om_user_color())
			butcolor['menu'].entryconfigure('true colors',command=lambda:set_om_color('true colors','#d8d8d8',1))
			set_om_color(*allcol[gfx.nmol%6]) #initial value
			butadd2lst = tk.Button(frame1, text='Add' ,width=10,height=1,command=lambda :mod.add_mol(gfx,self.mollist2,nameentry,self.modlist2,self.mollist,butcolor,rendvar,truecolor,varcolor,self.mollist3))
			butadd2lst.place( x=4,y=170)
			butrm2lst = tk.Button(frame1, text='Remove',width=10,height=1, command=lambda :mod.remove_mol(gfx,self.mollist2,self.mollist,self.mollist3))
			butrm2lst.place(x=230,y=170)
			butup = tk.Button(frame1, text='⋀\n⋀',width=1, command=lambda :mod.up_or_down_mol(gfx,self.mollist2,self.mollist,self.mollist3,1))
			butup.place(x=350,y=25)
			butdown = tk.Button(frame1, text='⋁\n⋁',width=1, command=lambda :mod.up_or_down_mol(gfx,self.mollist2,self.mollist,self.mollist3,-1))
			butdown.place(x=350,y=110)
			txtlb = tk.Label(frame1, text='Molecules List :')
			txtlb.place(x=230,y=1)
			self.mollist2=tk.Listbox(frame1,selectmode=tk.EXTENDED,width=12, height=8) #,background='black',foreground='white')
			scrollbar2 = tk.Scrollbar(frame1, orient=tk.VERTICAL)
			scrollbar2.config(command=self.mollist2.yview)
			self.mollist2.config(yscrollcommand=scrollbar2.set)
			scrollbar2.place(x= 331, y=26,width=14,height=129)
			self.mollist2.place(x=230,y=25)
			for mol in gfx.mol:
				self.mollist2.insert(tk.END,mol.un)
				self.mollist2.itemconfigure(tk.END,background=mod.vtk2tkhex_color(mol.col))
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
		#                        CONSTELLATION                             #
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
			frame2 = tk.Frame(nb)
			txtslb = tk.Label(frame2, text='Select Molecule :')
			txtslb.place( x=4,y=1)
			self.mollist3=tk.Listbox(frame2,selectmode=tk.EXTENDED,width=12, height=11)
			scrollbar = tk.Scrollbar(frame2, orient=tk.VERTICAL)
			scrollbar.config(command=self.mollist3.yview)
			self.mollist3.config(yscrollcommand=scrollbar.set)
			scrollbar.place(x= 105, y=26,width=14,height=178)
			self.mollist3.place(x=4,y=25)
			for mol in gfx.mol:
				self.mollist3.insert(tk.END,mol.un)
				self.mollist3.itemconfigure(tk.END,background=mod.vtk2tkhex_color(mol.col))
			self.mollist3.select_set(0,tk.END)
			def select_all_const_list():
				self.mollist3.select_set(0,tk.END)
			nb.bind('<<NotebookTabChanged>>', lambda evt: select_all_const_list())
			butselall = tk.Button(frame2, text='Select All', command=lambda :select_all_const_list())
			butselall.place( x=130,y=20)
			varindi = tk.StringVar(frame2)
			varcoli = tk.StringVar(frame2)
			if gfx.sym!=None:
				varindi.set(gfx.sym.indi)
			else :
				varindi.set(1)
			if gfx.sym!=None:
				varcoli.set(gfx.sym.coli)
			else :
				varcoli.set(1)
			butcoli= tk.Checkbutton(frame2, text="Collision Filter",variable=varcoli,onvalue=1)
			butcoli.place(x=135,y=60)
			cbig= tk.Checkbutton(frame2, text="Individual Constellation",variable=varindi,onvalue=1)
			cbig.place(x=135,y=85)
			txtsmn = tk.Label(frame2, text='Minimal Distance \n to Box Edges :')
			txtsmn.place( x=135,y=115)
			mdbeentry = tk.Entry(frame2, width=8,background='beige') #Minimal distance to box edges
			mdbeentry.place(  x=245,y=128)
			if gfx.sym!=None:
				mdbeentry.delete(0,tk.END)
				mdbeentry.insert(tk.END,str(gfx.sym.mdbe))
			else :
				mdbeentry.insert(tk.END,str(0))
			butbldconst = tk.Button(frame2, text='Build Constellation', command=lambda :load_sym(gfx))
			butbldconst.place( x=170,y=175)
			def load_sym(gfx):
				gfx.sym=sym.Sym(self.mollist3,gfx,float(mdbeentry.get()),varindi.get(),varcoli.get())
				gfx.sym.main()
			'''
			def bind_return(event):
				if nb.index("current")==0:
					butadd2lst.invoke()
				elif nb.index("current")==1:
					butbldconst.invoke()
			frame.bind("<Return>",bind_return)
			'''
			nb.add(frame1, text='Molecules', padding=1)
			nb.add(frame2, text='Constellation', padding=1)
			nb.pack(expand=1, fill='both')
			self.molwizopen=1
		else:
			MB.showwarning('Info','The molecule wizard is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.molwizquit(frame))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def molwizquit(self,frame):
		self.modlist2.unbind_all('<ButtonRelease>')
		self.modlist2=None
		self.molwizopen=0
		frame.destroy()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              FITTING SETUP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def setup_wiz(self,gfx):
		if not self.setupwizopen:
			frame = tk.Toplevel(width=250 ,height=120)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.title('Fitting Setup')
			frame.transient(self.root)
			#---
			reslab = tk.Label(frame, text='Resolution ')
			reslab.place( x=5,y=6)
			#---
			limitlab = tk.Label(frame, text='(limit =             ) :')
			limitlab.place( x=85,y=6)
			varlim=tk.StringVar()
			limlabel = tk.Label(frame,textvariable=varlim,width=6)
			limlabel.place( x=130,y=6)
			if gfx.map!=[]:
				(resmin, resmax) = gfx.map[0].map_res_range()
				varlim.set('%.2f'%resmax)
			else :
				varlim.set('NA')
			#----
			reslowlab = tk.Label(frame, text='Low :')
			reslowlab.place( x=32,y=30)
			reslowentry = tk.Entry(frame,width=6, background='beige')
			reslowentry.place(x=68,y=30)
			if gfx.fit!=None:
				reslowentry.delete(0,tk.END)
				reslowentry.insert(tk.END,str(gfx.fit.reslow))
			else :
				reslowentry.insert(tk.END,'')
			#----
			reshighlab = tk.Label(frame, text='High :')
			reshighlab.place( x=140,y=30)
			reshighentry = tk.Entry(frame,width=6, background='beige')
			reshighentry.place(x=180,y=30)
			if gfx.fit!=None:
				reshighentry.delete(0,tk.END)
				reshighentry.insert(tk.END,str(gfx.fit.reshigh))
			else :
				reshighentry.insert(tk.END,'')
			#---
			txtne = tk.Label(frame, text='EM Scale :')
			txtne.place( x=114,y=60)
			self.fitscaleentry = tk.Entry(frame, width=6,background='beige')
			self.fitscaleentry.place(  x=180,y=60)
			if gfx.map!=[]:
				self.fitscaleentry.delete(0,tk.END)
				self.fitscaleentry.insert(tk.END,str(gfx.map[0].scale))
			else :
				self.fitscaleentry.insert(tk.END,"1.0")
			#---
			vartest = tk.StringVar(frame)
			if gfx.fit!=None:
				vartest.set(gfx.fit.vartest)
			else :
				vartest.set(0)
			cbtest= tk.Checkbutton(frame, text="Test",variable=vartest,onvalue=1)
			cbtest.place(x=145,y=90)
			cbtest.place_forget()
			#plus de res low, definition par defaut :
			#reslow = 300
			butcfg = tk.Button(frame, text='Configure', command=lambda : uro.configure_fitting(gfx,reslowentry.get(),reshighentry.get(),vartest.get(),float(self.fitscaleentry.get())))
			butcfg.place(x=4,y=85)
			'''
			def bind_return(event):
				butcfg.invoke()
			frame.bind("<Return>",bind_return)
			'''
			self.setupwizopen=1
		else:
			MB.showwarning('Info','The Fitting Setup wizard is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.setupwizquit(frame))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def setupwizquit(self,frame):
		self.setupwizopen=0
		frame.destroy()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Propagate Scale
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def propagate_scale(self,gfx,scale,caller):
		if caller == 'fit':
			if self.mapwizopen:
				self.mapscaleentry.delete(0,tk.END)
				self.mapscaleentry.insert(tk.END,str(scale))
				self.sgv.set('%.2f'%(gfx.ps.radius * scale))
				self.dzentry.delete(0,tk.END)
				self.dzentry.insert(tk.END,'%.3f'%gfx.ps.dz) #pour dz : juste mise a jour pas besoin de scaler
				self.centry.delete(0,tk.END)
				self.centry.insert(tk.END,str(gfx.ps.c))
		if caller == 'map':
			if self.setupwizopen:
				self.fitscaleentry.delete(0,tk.END)
				self.fitscaleentry.insert(tk.END,str(scale))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              URO REFINEMENT FITING WIZ
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def fit_wiz(self,gfx):
		if not self.fitwizopen:
			frame = tk.Toplevel(width=210,height=250)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.title('URO Refinement')
			frame.transient(self.root)
			reslab = tk.Label(frame, text='Working Resolution Range :')
			reslab.place( x=5,y=5)
			#---
			reslowlab = tk.Label(frame, text='Low :')
			reslowlab.place( x=5,y=30)
			fitreslowentry = tk.Entry(frame,width=6, background='beige')
			fitreslowentry.place(x=40,y=30)
			if gfx.fit!=None:
				fitreslowentry.delete(0,tk.END)
				fitreslowentry.insert(tk.END,str(gfx.fit.fitreslow))
			else :
				fitreslowentry.insert(tk.END,str(400))
			#----
			reshighlab = tk.Label(frame, text='High :')
			reshighlab.place( x=110,y=30)
			fitreshighentry = tk.Entry(frame,width=6, background='beige')
			fitreshighentry.place(x=150,y=30)
			if gfx.fit!=None:
				fitreshighentry.delete(0,tk.END)
				fitreshighentry.insert(tk.END,str(gfx.fit.fitreshigh))
			else :
				fitreshighentry.insert(tk.END,str(20))
			#----
			nbclab = tk.Label(frame, text='Nb Cycles :')
			nbclab.place( x=50,y=55)
			nbcentry = tk.Entry(frame,width=6, background='beige')
			nbcentry.place(x=150,y=55)
			if gfx.fit!=None:
				nbcentry.delete(0,tk.END)
				nbcentry.insert(tk.END,str(gfx.fit.nbc))
			else :
				nbcentry.insert(tk.END,len(gfx.mol))
			#----
			nbitlab = tk.Label(frame, text='Nb Iterations :')
			nbitlab.place( x=50,y=80)
			nbitentry = tk.Entry(frame,width=6, background='beige')
			nbitentry.place(x=150,y=80)
			if gfx.fit!=None:
				nbitentry.delete(0,tk.END)
				nbitentry.insert(tk.END,str(gfx.fit.nbit))
			else :
				nbitentry.insert(tk.END,str(10))
			#----
			varbfact = tk.StringVar(frame)
			if gfx.fit!=None:
				varbfact.set(gfx.fit.varbfact)
			else :
				varbfact.set(0)
			cbbfact= tk.Checkbutton(frame, text="Use B Factor",variable=varbfact,onvalue=1)
			cbbfact.place(x=95,y=110)
			#----
			butrfhconst = tk.Button(frame,bg='red', text='Update Sym-Mates', command=lambda :refresh_all_const(gfx))
			butrfhconst.place(x=35,y=140)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def refresh_all_const(gfx):
				if gfx.sym!=None:
					gfx.sym.refresh_all_constelation(gfx)
				else :
					MB.showwarning('Info','No Constellation has been loaded')
			#----
			butfit = tk.Button(frame, text='Refine', command=lambda : uro.refinement(gfx,fitreslowentry.get(),fitreshighentry.get(),nbcentry.get(),nbitentry.get(),varbfact.get(),res))
			butfit.place(x=4,y=175)
			'''
			def bind_return(event):
				butfit.invoke()
			frame.bind("<Return>",bind_return)
			'''
			#----
			butundo = tk.Button(frame, text='Undo', command=lambda : undo_or_redo('G0'))
			butundo.place(x=83,y=175)
			butredo = tk.Button(frame, text='Redo', command=lambda : undo_or_redo('G1'))
			butredo.place(x=147,y=175)
			#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def undo_or_redo(gfile):
				if gfx.fit != None:
					gfx.fit.assign_pv_from_file(gfile,res)
			#----Result
			resultslab = tk.Label(frame, text='------------------- Results -------------------')
			resultslab.place( x=0,y=210)
			res=tk.StringVar()
			reslabel = tk.Label(frame,textvariable=res,width=26,background='beige')
			reslabel.place( x=0,y=230)
			self.fitwizopen=1
		else:
			MB.showwarning('Info','The Fitting wizard is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.fitwizquit(frame))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def fitwizquit(self,frame):
		self.fitwizopen=0
		frame.destroy()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              INTERACTIVE FITING WIZ
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def int_fit_wiz(self,gfx):
		if not self.ifitwizopen:
			frame = tk.Toplevel(width=210,height=170)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.title('Interactive Fitting')
			frame.transient(self.root)
			reslab = tk.Label(frame, text='Working Resolution Range :')
			reslab.place( x=5,y=5)
			#---
			reslowlab = tk.Label(frame, text='Low :')
			reslowlab.place( x=5,y=30)
			self.fitreslowentry = tk.Entry(frame,width=6, background='beige')
			self.fitreslowentry.place(x=40,y=30)
			if gfx.fit!=None:
				self.fitreslowentry.delete(0,tk.END)
				self.fitreslowentry.insert(tk.END,str(gfx.fit.fitreslow))
			else :
				self.fitreslowentry.insert(tk.END,str(400))
			#----
			reshighlab = tk.Label(frame, text='High :')
			reshighlab.place( x=110,y=30)
			self.fitreshighentry = tk.Entry(frame,width=6, background='beige')
			self.fitreshighentry.place(x=150,y=30)
			if gfx.fit!=None:
				self.fitreshighentry.delete(0,tk.END)
				self.fitreshighentry.insert(tk.END,str(gfx.fit.fitreshigh))
			else :
				self.fitreshighentry.insert(tk.END,str(20))
			#----
			butrfhconst = tk.Button(frame,bg='red', text='Update Sym-Mates', command=lambda :refresh_all_const(gfx))
			butrfhconst.place(x=35,y=65)
			if gfx.ifit != None:
				self.fitreshighentry.configure(state=tk.DISABLED)
				self.fitreslowentry.configure(state=tk.DISABLED)
			else :
				self.fitreshighentry.configure(state=tk.NORMAL)
				self.fitreslowentry.configure(state=tk.NORMAL)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def refresh_all_const(gfx):
				if gfx.sym!=None:
					gfx.sym.refresh_all_constelation(gfx)
				else :
					MB.showwarning('Info','No Constellation has been loaded')
			#----
			butfit = tk.Button(frame, text='Start/Stop', command=lambda : ifit_launcher())
			butfit.place(x=4,y=100)
			#----Result
			resultslab = tk.Label(frame, text='------------------- Results -------------------')
			resultslab.place( x=0,y=134)
			res=tk.StringVar()
			reslabel = tk.Label(frame,textvariable=res,width=26,background='beige')
			reslabel.place( x=0,y=150)
			#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def ifit_launcher():
				if gfx.mol!=[]:
					if gfx.fit!=None:
						if gfx.ifit == None:
							self.fitreshighentry.configure(state=tk.DISABLED)
							self.fitreslowentry.configure(state=tk.DISABLED)
						else :
							self.fitreshighentry.configure(state=tk.NORMAL)
							self.fitreslowentry.configure(state=tk.NORMAL)
						uro.interactive_fitting(gfx,self.fitreslowentry.get(),self.fitreshighentry.get(),res)
					else :
						MB.showwarning('Info','Configure Fitting|Setup')
				else :
					MB.showwarning('Info','Go to Assemblage|Molecules and load molecules. Then configure Fitting|Setup')
			self.ifitwizopen=1
		else:
			MB.showwarning('Info','The Fitting wizard is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.ifitwizquit(frame,gfx))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def ifitwizquit(self,frame,gfx):
		if gfx.ifit != None:
			gfx.ifit.terminate_process()
		self.ifitwizopen=0
		frame.destroy()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              CC Profiler
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def cc_prof_wiz(self,gfx):
		if not self.ccprofwizopen:
			frame = tk.Toplevel(width=345,height=170)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.title('CC profiler')
			frame.transient(self.root)
			#---
			smplab = tk.Label(frame, text='Angular Sampling Range :')
			smplab.place( x=5,y=5)
			#---
			minlab = tk.Label(frame, text='Min :')
			minlab.place( x=5,y=30)
			minentry = tk.Entry(frame,width=6, background='beige')
			minentry.place(x=40,y=30)
			minentry.insert(tk.END,str(-180))
			#----
			maxlab = tk.Label(frame, text='Max :')
			maxlab.place( x=110,y=30)
			maxentry = tk.Entry(frame,width=6, background='beige')
			maxentry.place(x=150,y=30)
			maxentry.insert(tk.END,str(180))
			#----
			steplab = tk.Label(frame, text='Step :')
			steplab.place( x=5,y=60)
			stepentry = tk.Entry(frame,width=6, background='beige')
			stepentry.place(x=50,y=60)
			stepentry.insert(tk.END,str(5))
			#--
			varaxe = tk.StringVar(frame)
			varaxe.set('X') #initial value
			axelab = tk.Label(frame, text='Rotation Around Axe :')
			axelab.place(x=5,y=90)
			axeom=	tk.OptionMenu(frame,varaxe,'X','Y', 'Z')
			axeom.place(x=150,y=85)
			#---
			txtlb = tk.Label(frame, text='Select Molecule :')
			txtlb.place(x=230,y=5)
			mollist=tk.Listbox(frame,selectmode=tk.SINGLE,width=12, height=8)
			scrollbar = tk.Scrollbar(frame, orient=tk.VERTICAL)
			scrollbar.config(command=mollist.yview)
			mollist.config(yscrollcommand=scrollbar.set)
			scrollbar.place(x= 331, y=26,width=14,height=129)
			mollist.place(x=230,y=25)
			for mol in gfx.mol:
				mollist.insert(tk.END,mol.un)
				mollist.itemconfigure(tk.END,background=mod.vtk2tkhex_color(mol.col))
			#---
			butsample = tk.Button(frame, text='Sample Profile', command=lambda : uro.cc_profile(gfx,mollist,float(minentry.get()),float(maxentry.get()),float(stepentry.get()),varaxe.get()))
			butsample.place(x=5,y=130)
			self.ccprofwizopen=1
		else:
			MB.showwarning('Info','The CC Profiler wizard is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.ccprofwizquit(frame))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def ccprofwizquit(self,frame):
		self.ccprofwizopen=0
		frame.destroy()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              NCS_RMS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def ncs_rms_wiz(self,gfx):
		if not self.ncsrmswizopen:
			frame = tk.Toplevel(width=350,height=180)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.title('RMS Threshold')
			frame.transient(self.root)
			#---
			txtlb = tk.Label(frame, text='Select Molecule :')
			txtlb.place(x=4,y=5)
			self.mollist4=tk.Listbox(frame,selectmode=tk.SINGLE,width=12, height=8)
			scrollbar = tk.Scrollbar(frame, orient=tk.VERTICAL)
			scrollbar.config(command=self.mollist4.yview)
			self.mollist4.config(yscrollcommand=scrollbar.set)
			scrollbar.place(x= 105, y=26,width=14,height=129)
			self.mollist4.place(x=4,y=25)
			for mol in gfx.mol:
				self.mollist4.insert(tk.END,mol.un)
				self.mollist4.itemconfigure(tk.END,background=mod.vtk2tkhex_color(mol.col))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			def display_max_rmsd(event):
				varlim.set(uro.get_max_rmsd(gfx,self.mollist4))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			self.mollist4.bind_all('<ButtonRelease>', display_max_rmsd)
			#---
			varlim=tk.StringVar()
			rmslab = tk.Label(frame, text='RMS (Max Rot             ) :')
			rmslab.place( x=130,y=20)
			limlabel = tk.Label(frame,textvariable=varlim,width=6,anchor=tk.E)
			limlabel.place( x=215,y=20)
			varlim.set('NA')
			#--
			rmsentry = tk.Entry(frame,width=6, background='beige')
			rmsentry.place(x=285,y=20)
			rmsentry.insert(tk.END,'')
			if gfx.ncs_rms!=None:
				rmsentry.delete(0,tk.END)
				rmsentry.insert(tk.END,str(gfx.ncs_rms.rms))
			#----
			shiftlab = tk.Label(frame, text='Shift Type :')
			shiftlab.place( x=130,y=115)
			varstype = tk.StringVar(frame)
			varstype.set('Rotation') #initial value
			stypeom=tk.OptionMenu(frame,varstype,'Rotation','Translation', 'R + T')
			stypeom.place(x=210,y=110)
			if gfx.ncs_rms!=None:
				if gfx.ncs_rms.stype=='R':
					varstype.set('Rotation')
				elif gfx.ncs_rms.stype=='T':
					varstype.set('Translation')
				else :
					varstype.set('R + T')
			#--
			nbtrialslab = tk.Label(frame, text='Nb Trials :')
			nbtrialslab.place(x=130,y=50)
			trialentry = tk.Entry(frame,width=6, background='beige')
			trialentry.place(x=220,y=50)
			trialentry.insert(tk.END,str(30))
			if gfx.ncs_rms!=None:
				trialentry.delete(0,tk.END)
				trialentry.insert(tk.END,str(gfx.ncs_rms.nbtrials))
			#--
			nbitelab = tk.Label(frame, text='Nb Iterations :')
			nbitelab.place(x=130,y=80)
			iteentry = tk.Entry(frame,width=6, background='beige')
			iteentry.place(x=220,y=80)
			iteentry.insert(tk.END,str(100))
			if gfx.ncs_rms!=None:
				iteentry.delete(0,tk.END)
				iteentry.insert(tk.END,str(gfx.ncs_rms.nbite))
			#---
			butstart = tk.Button(frame, text='Start', command=lambda : uro.launch_rms_threshold(gfx,rmsentry.get(),varstype.get(),trialentry.get(),iteentry.get(),self.mollist4))
			butstart.place(x=130,y=145)
			self.ncsrmswizopen=1
		else:
			MB.showwarning('Info','The RMS Threshold wizard is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.ncsrmswizquit(frame))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def ncsrmswizquit(self,frame):
		self.ncsrmswizopen=0
		self.mollist4.unbind_all('<ButtonRelease>')
		self.mollist4 = None
		frame.destroy()	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            MAGNIFICATION FINDER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def magnum_wiz(self,gfx):
		if not self.magnumwizopen:
			frame = tk.Toplevel(width=260 ,height=75)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.title('Magnification Finder')
			frame.transient(self.root)
			#---
			steplab = tk.Label(frame, text='Step :')
			steplab.place( x=5,y=10)
			#----
			stepentry = tk.Entry(frame,width=6, background='beige')
			stepentry.place(x=55,y=10)
			stepentry.delete(0,tk.END)
			stepentry.insert(tk.END,str(0.01))
			#----
			nbsteplab = tk.Label(frame, text='Nb Steps :')
			nbsteplab.place( x=120,y=10)
			nbstepentry = tk.Entry(frame,width=6, background='beige')
			nbstepentry.place(x=195,y=10)
			nbstepentry.delete(0,tk.END)
			nbstepentry.insert(tk.END,str(5))
			butstart = tk.Button(frame, text='Start',width=9, command=lambda : uro.launch_magnum(gfx,stepentry.get(),nbstepentry.get()))
			butstart.place(x=85,y=40)
			'''
			def bind_return(event):
				butstart.invoke()
			frame.bind("<Return>",bind_return)
			'''
			self.magnumwizopen=1
		else:
			MB.showwarning('Info','The Magnification Finder is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.magnumwizquit(frame))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def magnumwizquit(self,frame):
		self.magnumwizopen=0
		frame.destroy()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Normal Modes Wizard (NMW)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def normal_modes_wiz(self,gfx):
		if not self.nmwizopen:
			frame = tk.Toplevel(width=388,height=235)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.geometry("388x235+600+120")
			frame.title('Normal Modes Wizard')
			frame.transient(self.root)
			nb = ttk.Notebook(frame, width=270, height=180)
			nb.pressed_index = None
			frame1=tk.Frame(nb)
			txtsmn = tk.Label(frame1, text='Select Model :')
			txtsmn.place( x=4,y=2)
			self.modlist3=tk.Listbox(frame1,selectmode=tk.SINGLE ,width=12, height=11)
			scrollbar = tk.Scrollbar(frame1, orient=tk.VERTICAL)
			scrollbar.config(command=self.modlist3.yview)
			self.modlist3.config(yscrollcommand=scrollbar.set)
			scrollbar.place(x= 105, y=26,width=14,height=178)
			self.modlist3.place(x=4,y=25)
			for mod in gfx.mod:
				self.modlist3.insert(tk.END,mod.un)
			butcptnm = tk.Button(frame1, text='Compute NM ->' ,width=10,height=1,command=lambda :Nma.compute_nm(gfx,self.modlist3,self.nmmodlist,self.nmmollist))
			butcptnm.place( x=130,y=60)
			butloadnm = tk.Button(frame1, text='Load NM File ->' ,width=10,height=1,command=lambda :Nma.load_nm(gfx,self.modlist3,self.nmmodlist,self.nmmollist))
			butloadnm.place( x=130,y=90)
			butview = tk.Button(frame1, text='View',width=3,height=1, command=lambda : viewer_laucher(gfx))
			butview.place(x=260,y=146)
			butsave = tk.Button(frame1, text='Save',width=3,height=1, command=lambda : Nma.save_nm(gfx,self.nmmodlist))
			butsave.place(x=315,y=146)
			butdel = tk.Button(frame1, text='Remove',width=10,height=1, command=lambda : Nma.remove_mod(gfx,self.nmmodlist,self.nmmollist))
			butdel.place(x=260,y=178)
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
			def viewer_laucher(gfx):
				for modo in gfx.mod:
					if self.nmmodlist.get(tk.ACTIVE) == modo.un:
						if not Nma.matrix_is_ready(gfx,modo):
							MB.showwarning('Info','Normal Modes are still computing ... please wait')
							return
						self.nm_viewer(gfx,modo)
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
			txtlb = tk.Label(frame1, text='NM Modlist :')
			txtlb.place(x=260,y=2)
			self.nmmodlist=tk.Listbox(frame1,selectmode=tk.EXTENDED,width=12, height=7)
			scrollbar2 = tk.Scrollbar(frame1, orient=tk.VERTICAL)
			scrollbar2.config(command=self.nmmodlist.yview)
			self.nmmodlist.config(yscrollcommand=scrollbar2.set)
			scrollbar2.place(x= 361, y=26,width=14,height=113)
			self.nmmodlist.place(x=260,y=25)
				
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
		#                            ANALYSES                              #
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
			frame2 = tk.Frame(nb)
			txtslb = tk.Label(frame2, text='Select Molecule :')
			txtslb.place( x=4,y=2)
			self.nmmollist=tk.Listbox(frame2,selectmode=tk.EXTENDED,width=12, height=11)
			scrollbar = tk.Scrollbar(frame2, orient=tk.VERTICAL)
			scrollbar.config(command=self.nmmollist.yview)
			self.nmmollist.config(yscrollcommand=scrollbar.set)
			scrollbar.place(x= 105, y=26,width=14,height=178)
			self.nmmollist.place(x=4,y=25)
			
			modenbtxt = tk.Label(frame2, text='Mode Number :')
			modenbtxt.place(x=135,y=50)
			modenbentry = tk.Entry(frame2, width=5,background='beige')
			modenbentry.place(x=260,y=50)
			modenbentry.delete(0,tk.END)
			modenbentry.insert(tk.END,str(7))
			
			steptxt = tk.Label(frame2, text='Step (rms/nat) :')
			steptxt.place(x=135,y=75)
			stepentry = tk.Entry(frame2, width=5,background='beige')
			stepentry.place(x=260,y=75)
			stepentry.delete(0,tk.END)
			stepentry.insert(tk.END,str(0.15))
			
			nbsteptxt = tk.Label(frame2, text='Number of Steps :')
			nbsteptxt.place(x=135,y=100)
			nbstepentry = tk.Entry(frame2, width=5,background='beige')
			nbstepentry.place(x=260,y=100)
			nbstepentry.delete(0,tk.END)
			nbstepentry.insert(tk.END,str(25))
			
			self.replacebut = tk.Button(frame2, text='Replace in Situ',width=15,height=1, command=lambda :in_situ_laucher())
			self.replacebut.place( x=170,y=135)
			self.acceptbut = tk.Button(frame2, text='Acceptance',width=15,height=1, command=lambda : self.nm_accept(gfx))
			self.acceptbut.place( x=170,y=165)
			self.acceptbut.configure(state=tk.DISABLED)
			def in_situ_laucher():
				for molo in gfx.mol:
					if self.nmmollist.get(tk.ACTIVE) == molo.un:
						Nma.replace_in_situ(gfx,molo,modenbentry.get(),stepentry.get(),nbstepentry.get())
						self.acceptbut.configure(state=tk.NORMAL)
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
			nb.add(frame1, text='Setup', padding=1)
			nb.add(frame2, text='Analysis', padding=1) #TXT 
			nb.pack(expand=1, fill='both')
			self.nmwizopen=1
			Nma.build_lists(gfx,self.nmmodlist,self.nmmollist)#apres open puisque test dans function sur if open
		else:
			MB.showwarning('Info','The Normal Modes Wizard is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.nmwizquit(frame,gfx))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def nmwizquit(self,frame,gfx):
		if not self.nmacceptopen:
			self.modlist3=None
			self.nmwizopen=0
			if gfx.nma != None:
				gfx.nma.stop_insitu_loop(gfx)
				gfx.nma.reset_normal_state(gfx)
			frame.destroy()
		else :
			MB.showwarning('Info','Close Acceptance Wizzard before Normal Modes Wizard')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  ACCEPTANCE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def nm_accept(self,gfx):
		if not self.nmacceptopen:
			frame = tk.Toplevel(width=500,height=250)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.geometry("+700+120")
			frame.transient(self.root)
			frame.title('Acceptance Wizard')
			gfx.nma.stop_insitu_loop(gfx)
			self.varconformerSitu.set(gfx.nma.oldon-gfx.nma.nbstep)
			stepnblb = tk.Label(frame, text='Step|Conformer Number')
			stepnblb.place( x=170,y=4)
			self.confscrollSitu = tk.Scale(frame,command=lambda o,g=gfx : g.nma.display_conformers_manual(o),variable=self.varconformerSitu,from_=-gfx.nma.nbstep,to=gfx.nma.nbstep,length=480,showvalue='yes',orient='horizontal')
			self.confscrollSitu.place(x=10,y=20)
			tiret1lb = tk.Label(frame, text='___________________________________________________________________________________')
			tiret1lb.place( x=0,y=60)
			aslb = tk.Label(frame, text='Amplitude Search')
			aslb.place( x=190,y=90)
			itvlb = tk.Label(frame, text='Interval :')
			itvlb.place( x=4,y=120)
			#----
			lowlab = tk.Label(frame, text='Low :')
			lowlab.place( x=82,y=120)
			lowentry = tk.Entry(frame,width=6, background='beige')
			lowentry.place(x=118,y=120)
			lowentry.delete(0,tk.END)
			lowentry.insert(tk.END,str(-gfx.nma.nbstep))
			#----
			highlab = tk.Label(frame, text='High :')
			highlab.place( x=190,y=120)
			highentry = tk.Entry(frame,width=6, background='beige')
			highentry.place(x=230,y=120)
			highentry.delete(0,tk.END)
			highentry.insert(tk.END,str(gfx.nma.nbstep))
			#----
			self.startbut = tk.Button(frame, text='Start',width=12,height=1, command=lambda : gfx.nma.norma(gfx,lowentry.get(),highentry.get()))
			self.startbut.place( x=190,y=150)
			tiret2lb = tk.Label(frame, text='___________________________________________________________________________________')
			tiret2lb.place( x=0,y=180)
			butsave = tk.Button(frame, text='Save Conformer Number',height=1, command=lambda :gfx.nma.save_conf_as(gfx))
			butsave.place( x=160,y=210)
			cpltnfv=nfv=tk.StringVar()
			stepnblb = tk.Label(frame, textvariable=self.varconformerSitu,background='beige')
			stepnblb.place( x=350,y=215)
			self.nmacceptopen=1
		else :
			MB.showwarning('Info','The Acceptance Wizard is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.nmacceptquit(frame,gfx))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def nmacceptquit(self,frame,gfx):
		self.nmacceptopen=0
		self.acceptbut.configure(state=tk.DISABLED)
		gfx.nma.reset_normal_state(gfx)
		frame.destroy()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Normal Modes Viewer (NMW)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def nm_viewer(self,gfx,modo):
		if not self.nmvieweropen:
			frame = tk.Toplevel(width=500,height=625)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.geometry("+700+120")
			frame.transient(self.root)
			frame.title('Normal Modes Viewer')
			varnbmode = tk.IntVar(frame)
			varnbmode.set(7)
			nbmodelab = tk.Label(frame, text='Mode Number :')
			nbmodelab.place(x=5,y=6)
			nbmodeom=tk.OptionMenu(frame,varnbmode,'7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26')
			nbmodeom.place(x=110,y=4)
			steptxt = tk.Label(frame, text='Step (rms/nat) :')
			steptxt.place(x=172,y=6)
			stepentry = tk.Entry(frame, width=5,background='beige')
			stepentry.place(x=277,y=6)
			stepentry.delete(0,tk.END)
			stepentry.insert(tk.END,str(0.15))
			
			nbsteptxt = tk.Label(frame, text='Number of Steps :')
			nbsteptxt.place(x=330,y=6)
			nbstepentry = tk.Entry(frame, width=5,background='beige')
			nbstepentry.place(x=450,y=6)
			nbstepentry.delete(0,tk.END)
			nbstepentry.insert(tk.END,str(25))
			
			butdisp = tk.Button(frame, text='Display',width=10,height=1, command=lambda :Nma.display_nm(gfx,modo,varnbmode.get(),stepentry.get(),nbstepentry.get()))
			butdisp.place( x=200,y=45)
			self.varloop.set(0)
			butloop= tk.Checkbutton(frame, text="Loop",variable=self.varloop,onvalue=1,command=lambda g=gfx:g.nma.viewer.loop_thread())
			butloop.place(x= 350, y=50)
			gfx.nma.viewer = Nma.nmviewer(gfx,frame)
			self.varconformer.set(0)
			self.confscroll = tk.Scale(frame,command=lambda o,g=gfx : g.nma.viewer.display_conformers_manual(o),variable=self.varconformer,from_=0,to=0,length=480,showvalue='yes',orient='horizontal')
			self.confscroll.place(x=10,y=75)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			
			frame.geometry("500x625+700+120")
			self.nmvieweropen=1
		else :
			MB.showwarning('Info','The Normal Modes Viewer is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.nmviewerquit(frame,gfx))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def nmviewerquit(self,frame,gfx):
		self.nmvieweropen=0
		del gfx.nma.viewer.renwin
		del gfx.nma.viewer
		frame.destroy()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              mols bookkeeping
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def bookkeeping(self,gfx):
		if not self.bookkeepingopen:
			frame = tk.Toplevel(width=150,height=100)
			frame.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame.title('BookKeeping')
			frame.transient(self.root)
			mod.build_bookkeeping(gfx)
			st = ScrolledText(frame,width=200,height=20,background='white')
			st.pack (side='top', fill='both', expand=1)
			st.delete('0.0',tk.END)
			st.insert(tk.END,open(gfx.tmpdir + '/mols.bk','r').read())
			self.bookkeepingopen=1
		else:
			MB.showwarning('Info','BookKeeping is already open')
			return
		frame.protocol("WM_DELETE_WINDOW", lambda:self.bookkeepingquit(frame))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def bookkeepingquit(self,frame):
		self.bookkeepingopen=0
		frame.destroy()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              VISUAL CROP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def crop_wiz(self,gfx):#même interface que dans le wiz map sans le bouton crop map
		if not self.viscropopen :
			frame3 = tk.Toplevel(width=240,height=125)
			frame3.iconbitmap('@' + gfx.vedadir + 'doc/veda.xbm')
			frame3.title('Visual Crop')
			frame3.transient(self.root)
			butshow3 = tk.Button(frame3, text='Start', command=lambda : gfx.crop.show())
			butshow3.place(  x=50,y=5)
			butstop = tk.Button(frame3, text='Stop', command=lambda : gfx.crop.hide())
			butstop.place(  x=130,y=5)
			xllab = tk.Label(frame3, text='x low :')
			xllab.place( x=5,y=45)
			xlentry = tk.Entry(frame3,width=7, background='beige')
			xlentry.place(x=50,y=45)
			yllab = tk.Label(frame3, text='y low :')
			yllab.place( x=5,y=70)
			ylentry = tk.Entry(frame3,width=7, background='beige')
			ylentry.place(x=50,y=70)
			zllab = tk.Label(frame3, text='z low :')
			zllab.place( x=5,y=95)
			zlentry = tk.Entry(frame3,width=7, background='beige')
			zlentry.place(x=50,y=95)
			xulab = tk.Label(frame3, text='x up :')
			xulab.place( x=130,y=45)
			xuentry = tk.Entry(frame3,width=7, background='beige')
			xuentry.place(x=170,y=45)
			yulab = tk.Label(frame3, text='y up :')
			yulab.place( x=130,y=70)
			yuentry = tk.Entry(frame3,width=7, background='beige')
			yuentry.place(x=170,y=70)
			zulab = tk.Label(frame3, text='z up :')
			zulab.place( x=130,y=95)
			zuentry = tk.Entry(frame3,width=7, background='beige')
			zuentry.place(x=170,y=95)
			cropentry=(xlentry,xuentry,ylentry,yuentry,zlentry,zuentry)
			if gfx.crop!=None:
				xlentry.delete(0,tk.END)
				xlentry.insert(tk.END,str(gfx.crop.entryval[0]))
				xuentry.delete(0,tk.END)
				xuentry.insert(tk.END,str(gfx.crop.entryval[1]))
				ylentry.delete(0,tk.END)
				ylentry.insert(tk.END,str(gfx.crop.entryval[2]))
				yuentry.delete(0,tk.END)
				yuentry.insert(tk.END,str(gfx.crop.entryval[3]))
				zlentry.delete(0,tk.END)
				zlentry.insert(tk.END,str(gfx.crop.entryval[4]))
				zuentry.delete(0,tk.END)
				zuentry.insert(tk.END,str(gfx.crop.entryval[5]))
				gfx.crop=Map.Crop(gfx,gfx.map[0].iso,cropentry,gfx.crop.entryval)
		else :
			MB.showwarning('Info','Visual Crop is already open')
			return
		frame3.protocol("WM_DELETE_WINDOW", lambda:self.visualcropquit(frame3,gfx))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def visualcropquit(self,frame,gfx):
		self.viscropopen=0
		frame.destroy()
		if gfx.crop != None:
			gfx.crop.hide()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################################END##########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
