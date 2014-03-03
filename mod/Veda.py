# -*- coding: utf-8 -*-
# Main Module of VEDA
#
# Copyright:
# 2009-01-01 Gael Goret  gael.goret@ibs.fr
#
# Last modified:
# 2010-01-25 Gael Goret  gael.goret@ibs.fr
try :
	import Tkinter as tk
except :
	import tkinter as tk
import vtk
import gfx,itf
root = tk.Tk()
root.geometry("+0+0")
frame=itf.create_frame(root)
gfx=gfx.gfx()
gfx.setenv(root)
itfenv=itf.itf(root,gfx)
itfenv.init_menu(root,gfx)
itfenv.init_mod_list(frame,gfx)
root.mainloop()
