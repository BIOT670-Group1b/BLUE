import tkinter as tk
import samfile

from tkinter import *

counter = 0

class InitDNAFiles:

    def __init__(self, filepath, type):
        self.filepath = filepath
        self.type = type

    def init_chr_file(self):
        try:
            self.type = self.type.replace("CHR", "")  # To access the FASTA, we need the exact chromosome name
            samfile.load_chromosome(self.filepath,self.type)
        except:
            raise Exception("Error loading chromosome")
        finally:
            root.quit()
            root.destroy()

    def init_sam_file(self):
        try:
            loading_msg("\n\nLoading BAM file...")
            samfile.load_sam(self.filepath)
            loading_msg("\n\nStoring BAM data...")
            samfile.store_all_reads()
            loading_msg("\n\nLocating fusions...")
            samfile.map_fusion_reads()
        except:
            raise Exception("Error loading bam file")
        finally:
            root.quit()
            root.destroy()

    def init_gtf_file(self):
        try:
            samfile.load_genome(self.filepath)
        except:
            raise Exception("Error loading the gtf file")
        finally:
            root.quit()
            root.destroy()

    def create_loading(self):
        global root
        root = tk.Tk()
        root.title("Loading")
        root.geometry("250x100")

        if "CHR" in self.type:
            label = tk.Label(root, text="\n\nLoading chromosome...\n")
            label.pack()
            root.after(200, self.init_chr_file)
        if self.type == "SAM":
            label = tk.Label(root, text="\n\nLoading BAM file...\n")
            label.pack()
            root.after(200, self.init_sam_file)
        if self.type == "GTF":
            label = tk.Label(root, text="\n\nLoading/creating reference database..")
            label.pack()
            root.after(200, self.init_gtf_file)

        root.resizable(False,False)
        root.mainloop()

def loading_msg(msg):
    var1 = StringVar()
    label = tk.Label(root, textvariable=var1)
    label.pack()
    var1.set(msg)