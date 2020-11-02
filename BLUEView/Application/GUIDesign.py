import samfile
import tkinter as tk
import InitDNAFiles
import webbrowser
import threading

from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox

# GUI related funtions
def on_release(event):
    event.widget.config(cursor="")
def on_enter(event):
    event.widget.config(cursor="")

# Move Fusion/Chromosome viewers with mouse
# Chromosomal viewer gets position and adjusts due to the large sequence available
def move_start(event): # Get the current position upon clicking in the viewer
    event.widget.scan_mark(event.x, event.y)
    event.widget.config(cursor="hand1")
def move_move(event): # Adjust the chromosome sequence when nearing the current portions beginning/end
    chrom_indx = chrom_area.index(INSERT) # Default 1.0 if no chromosome is loaded
    line_num = int(str(chrom_indx)[0:1]) # Get the line number
    # Get the position in the viewer to update if reaches above or below a certain threshold
    viewer_pos = int(str(chrom_indx)[2:])
    chrom_seq = chrom_area.get(chrom_indx, END)
    chrom_seq = chrom_seq.strip(" \n")

    try: CHROMOSOME_SEQ
    except NameError: pass
    else: relative_pos = CHROMOSOME_SEQ.find(chrom_seq) # Find the chromosome position only if initialized

    event.widget.scan_dragto(event.x, 0) # This moves the viewer by dragging it

    if line_num==3: # Focus on th line with the sequence
        if viewer_pos<=1000 or viewer_pos>=9000: # Check if the viewer is reaching threshold of the loaded sequence
            # Given the current sequence, find where it is relative to the chromosome
            view_diff = relative_pos - viewer_pos
            if view_diff > 0: # If upon reaching the threshold we still have more bases to view, load them up
                if relative_pos >= VIEWER_SPAN:
                    chrom_area.unbind_class('post-class-bindings', '<ButtonPress-1>')
                    chrom_area.unbind_class('post-class-bindings', '<B1-Motion>')  # Deactive so sequence can properly update
                    win.after(500, rebind)  # Delay rebinding to allow time for viewer to repopulate
                    new_start_span = int(relative_pos - (VIEWER_SPAN/2)) + 50 # Adjust for position to remain in the same place
                    new_end_span = int(relative_pos + (VIEWER_SPAN/2)) + 50

                    reset_chromosome() # Reset the chromosome area and populate with the new section
                    insert_chrom_seq(CHROMOSOME_SEQ, new_start_span, new_end_span, True)

# Reactivates unbound functions
def rebind():
    global bind_chrom_move
    global bind_chrom_start
    bind_chrom_start = chrom_area.bind_class('post-class-bindings', '<ButtonPress-1>', move_start)
    bind_chrom_move = chrom_area.bind_class('post-class-bindings','<B1-Motion>', move_move)

# Displays the genes currently in view of the chromosome viewer
def update_gene_area():
    try: SEL_CHROMOSOME # First check if a chromosome has been loaded
    except NameError: pass
    else:
        try: genome_open # Then check if a reference genome has been loaded (for obtaining names)
        except NameError: pass
        else:
            genes_in_area.delete("1.0", END) # Clear the area
            chrom_indx = chrom_area.index("@0,0") # Get the first visible position in the viewer
            chrom_indx = chrom_indx.replace("1.","3.") #Adjust the lines to accurately represent the position
            chrom_indx = chrom_indx.replace("2.", "3.")
            chrom_seq = chrom_area.get(chrom_indx, END) # Get the sequence
            chrom_seq = chrom_seq.strip(" \n")
            relative_pos = CHROMOSOME_SEQ.find(chrom_seq) # Find the position the viewer sequence is relative to the chromosome
            genes = samfile.get_gene_info(relative_pos, SEL_CHROMOSOME) # Use the start and chromosome to get the info
            genes_in_area.config(state=NORMAL)
            genes_in_area.delete("1.0", END)
            for gene in genes:
                gene_string = gene.gene_name + ": " + str(gene.start) + " - " + str(gene.end)
                genes_in_area.insert(INSERT, gene_string + "\n")
            genes_in_area.tag_add("center", "1.0", END)
            genes_in_area.config(state=DISABLED)

    win.after(2000, update_gene_area) # As long as the chromosome and genome is loaded, search for genes every 2 seconds

# Go to the position or gene position entered in the search bar
def search_position(entered):
    entered = entered.replace(",","") # Format the position in case of commas

    try:
        CHROMOSOME_SEQ # Check first to see if the chromosome has been loaded
    except NameError:
        pass
    else:
        try: entered = int(entered) # Check to see if we need to find the sequence based on name or position alone
        except:
            try: genome_open # If the reference is loaded, then we can locate position by gene name
            except: pass
            else:
                if isinstance(entered, str) and entered != '': # If a string and not empty, locate the potential gene
                    entered = samfile.get_start_position_by_name(entered)
                    entered = int(entered) + 50 # Adjust for half length of viewer

        if not isinstance(entered, str) and entered != '': # If we get to here with a numeric value, go to the position
            chrom_region_start = int(entered) - 5000 # e.g. 75,758,257 - 5000 = 75,753,257
            chrom_region_end = int(entered) + 5000
            chrom_area.config(state=NORMAL)
            chrom_area.delete("1.0", END)  # Reset the viewer with the new positioned reads
            try: insert_chrom_seq(CHROMOSOME_SEQ, chrom_region_start, chrom_region_end, True) # No error if out of bounds, just blank
            except: pass
        elif entered != '': display_msg("Please load the reference file before searching by gene name.")

    try:
        sam_open # Check if the sam file has been loaded
    except NameError:
        pass
    else: # If so, go to the position in the fusion viewer if exists
        fuse_pos, fuse_len = samfile.get_fusions()
        try: adjusted_search_pos = int(entered) - int(fuse_pos[0]) # No error thrown if position is out of bounds
        except: pass
        else: fusion_area.see("1." + str(adjusted_search_pos+8))


# Looks for text and highlights it
def search(text_widget, keyword, tag):
    pos = '1.0'
    while True:
        idx = text_widget.search(keyword,pos,END)
        if not idx:
            break
        pos = '{}+{}c'.format(idx,len(keyword))
        text_widget.tag_add(tag,idx,pos)

# Display the document that shows how to use the application and what everything means
def open_how_to():
    WinUse = Tk()
    WinUse.title("How to Use Application")

    frameH = Frame(WinUse)
    frameH.pack()

    helpBox = Text(frameH, width=70, height=29, wrap=WORD)
    helpBox.grid(row=0)

    HelpText1 = '''BAM ONLY: To use this application, simply select a BAM file from the 
                  'File' dropdown and the program will read it. If a fusion was found, 
                  the reads will be displayed in the Fusion Area in the middle of the 
                  application. The fusion gene will be displayed and color coded on the 
                  left side, while it's supplement will show up on the right. Positions 
                  mark the start of the base that follows and each 'marker' is 5 bases 
                  in length.'''

    HelpText = "".join(HelpText1.splitlines()) + "\n\n"

    HelpText2 = '''BAM W/ ANNOTATIONS: In order for the program to retrieve the gene 
                  names and other relevant information, a GTF reference file will need 
                  to be loaded into the program. This can be done by selecting the 
                  reference the 'Ref' tab, the type of genome, and the file associated 
                  with that genome. The program comes preloaded with Ensembles GChr38 
                  annotation file for use. You can load the reference file before 
                  loading the BAM file and vice versa. For fusions found, the Fusion 
                  gene list and the Supplement gene list will populate with discovered 
                  fusions. More information about each gene can be found by clicking 
                  'More Information' under each gene name. The summary section will 
                  show key information about the fusion that was found. 
                  ***NOTE: The reference file will have to be reloaded if loading bam 
                  files one after another to populate the gene lists and summary area.'''

    HelpText = HelpText + "".join(HelpText2.splitlines())

    HelpText3 = '''CHROMOSOMES: Each chromosome (hg38) can be loaded by using the
                  dropdown located near the top of the application. Selecting a sequence
                  will load it into the Chrom-view. Once there, the user may drag left 
                  or right and genes will populate based on what is currently in view.
                  Users may also use the search function at the top to focus on a 
                  position or go to where a specific gene is located at.'''

    HelpText = HelpText + "\n\n" +  str("".join(HelpText3.splitlines()))
    helpBox.insert(INSERT, HelpText)

    WinUse.resizable(False,False)
    WinUse.mainloop()

# Send an email to the administrator to report any bugs/areas of improvement
def report_error():
    print("Nothing to see here.")

# Changes the color scheme of the application
def change_theme(theme):
    global DEFAULT_THEME
    DEFAULT_THEME = theme

    for wid in widget_list:
        try: wid.configure(bg=theme)
        except:
            try:
                wid.config(bg=theme)
            except: pass


# Given the zipped chromosome file, load the chromosome fasta into the chromosomal window
def load_chromosome(event):
    global CHROM_START # Starting position of the chromosome region. Default = 100000
    global CHROM_END # Ending position of the chromosome region. Default = 110000
    global VIEWER_SPAN # Number of bases that span the chromosomal viewer.
    global CHROMOSOME_SEQ
    global SEL_CHROMOSOME # Selected chromosome

    # Make sure these distances are 10,000 apart unless search_position distance numbers are other than 5000 (double for viewer span)
    CHROM_START = 10000
    CHROM_END = 20000
    VIEWER_SPAN = CHROM_END - CHROM_START

    SEL_CHROMOSOME = event.widget.get().replace('chr', '')
    file = "Chromosomes/Homo_sapiens.GRCh38.dna.chromosome." + str(SEL_CHROMOSOME) + ".fa.gz"
    if "Select" not in SEL_CHROMOSOME:
        try:
            reset_chromosome()
            type = "CHR" + str(SEL_CHROMOSOME) # Set the chromosome name for extraction
            init= InitDNAFiles.InitDNAFiles(file, type)
            init.create_loading()
            CHROMOSOME_SEQ = samfile.get_chromosome_seq()
            insert_chrom_seq(CHROMOSOME_SEQ, CHROM_START, CHROM_END, True)
        except Exception:
            display_error("Error loading chromosome. Please make sure you have the correct chromosome file: " + file)


# Open and Load BAM File
def _open_sam_file():
    global sam_open
    file = filedialog.askopenfilename(filetypes=(("bam files","*.bam"),("All files","*.*")))
    if file:
        try:
            reset_reads()  # Blank out the fusion area
            init = InitDNAFiles.InitDNAFiles(file, "SAM")
            init.create_loading() # Create a loading bar and load the bam/sam file
            sam_open = True
        except:
            display_error("Error loading BAM file. Please make sure the bai file exists in the same location as the BAM.")


# Open and Load Genome
def _open_genome_file(refGen):
    global fusion_genes
    global genome_open

    fusion_genes = []

    # Ask the user to manually select a reference file
    #file = filedialog.askopenfilename(filetypes=(("gtf files", "*.gtf"), ("All files", "*.*")))
    file = "Genomes/Homo_sapiens.GRCh38.100.gtf"
    if file:
        try:
            init = InitDNAFiles.InitDNAFiles(file, "GTF")
            init.create_loading() # Create the loading bar and load genome
            refLabel = Label(frameRef, text='Reference Genome ' + str(refGen) + " in Use", bg=DEFAULT_THEME)
            refLabel.grid(row=0, column=0)
            widget_list.append(refLabel) # Thematic changes list
            genome_open = True
            try: sam_open # Check to see if a sam file has been loaded, if so reset fusion view
            except NameError: pass
            else:
                reset_reads()
                samfile.map_fusion_reads()
        except:
            display_error("Error loading GTF file. The file may be downloaded at this address: "
                          "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz")


# Notify of error
def display_error(errMsg=""):
    messagebox.showerror("Error", errMsg)
# Display regular message
def display_msg(msg):
    messagebox.showinfo("Notice", msg)

# Clears the fusion gene area when a new file is loaded
def reset_reads():
    fusion_area.config(state=NORMAL)
    Gene1Txt.config(state=NORMAL)
    Gene2Txt.config(state=NORMAL)
    SummaryTxt.config(state=NORMAL)

    Gene1Txt.delete('1.0', END)
    Gene2Txt.delete('1.0', END)
    SummaryTxt.delete('1.0', END)
    try:
        Gene1Label
    except NameError:
        pass
    else:
        Gene1Label.grid_forget()
        Gene2Label.grid_forget()
    fusion_area.delete('1.0', END)

# Reset the chromosome area
def reset_chromosome():
    chrom_area.config(state=NORMAL)
    chrom_area.delete('1.0', END)
    chrom_area.config(state=DISABLED)


# Insert the chromosome sequence into the viewer at the bottom given a seq, start and end position
def insert_chrom_seq(seq, start, end, isMiddle):
    middle_indx = "3." + str(int(VIEWER_SPAN / 2))
    if not isMiddle: middle_indx = "3.0"

    # The chromosome view has fixed markers every 50 bases
    # Such we must adjust the sequence to display accurately if not divisible by 50
    start = round(start/50) * 50
    remainder = (start % 50) - 1
    seq = seq[start+remainder:end+remainder]
    seq_len = len(seq)

    chrom_area.config(state=NORMAL)
    positions = samfile.get_position_headers(seq_len,start, 50)
    create_headers(seq, positions, chrom_area)

    # Initialize the chromosome region sequence
    chrom_area.insert(INSERT, "\t" + seq)
    chrom_area.see(middle_indx) # Default to the middle to allow space for left and right view
    chrom_area.config(state=DISABLED)


# Fusion reads are inserted line by line into the fusion area
# Other views are created depending on loading of GTF file
def insert_read(fusion, supplement, fusion_name, fusion_id, supp_name, supp_id, isFirst,
                position_headers, fusion_position, supp_position, fusion_length, supp_length):
    fusion_area.config(state=NORMAL)
    # Create the headers if at least one fusion was found
    if isFirst and fusion != "No Fusions Detected":
        create_headers(fusion, position_headers, fusion_area)

    if fusion == "No Fusions Detected":
        fusion = "\n\n" + fusion
        fusion_area.tag_config('NA', foreground="indian red", font=("Georgia", "18", "bold"), justify='center')
        fusion_area.insert(INSERT, "\t" + fusion, 'NA')

    if supplement != "No Fusions Detected":
        fusion_area.insert(INSERT, "\t" + fusion)
        fusion_area.insert(INSERT, '\n')

        len_diff = fusion_length - supp_length
        fusion_area.tag_config(fusion, background='pale green')
        fusion_area.tag_config(supplement, background='light salmon')
        # Highlights the fusion and supplement gene
        search(fusion_area, fusion, fusion)
        search(fusion_area, supplement, supplement)

        focus_fusion_area(len_diff)

    # Lock the text box
    fusion_area.config(state=DISABLED)
    if fusion_name != "" and fusion_name not in fusion_genes:
        add_gene_info(fusion_name, supp_name, fusion_id, supp_id)
        add_fusion_summary(fusion_name, supp_name, fusion_position, supp_position, fusion_length, supp_length)
        fusion_genes.append(fusion_name)
        fusion_genes.append(supp_name)


# Place the position markers for the first row read
def create_headers(seq, positions, viewer):
    length = len(seq)
    dist_between = 9 # Number of markers before a positional marker is set
    pos_marker = "\u2502"
    mini_markers = ("    \u2758") * dist_between
    viewer.tag_config('Header', background="gray81")
    viewer.insert(INSERT, "|", "Marker") # Small font will shift the position markers slightly so to land in between bases
    viewer.insert(INSERT, " " * length, 'Header') # Initialize the header region so indeces may be referenced
    viewer.insert(INSERT, "\n")
    viewer.insert(INSERT, "|", "Marker") # Twice for the positions and the line markers themselves
    viewer.insert(INSERT, " " * length, 'Header')
    viewer.insert(INSERT, "\n")

    start_position = positions[0]

    for x in range(len(positions)):
        position_chars = len(str(positions[x])) # Get the number of characters in the position (as string)
        pos_adjusted_indx = round(position_chars/2)
        adjusted_start_indx = 8 - pos_adjusted_indx # Used to shift the position into the center of the marker (8 = length of tab)

        if x == 0:
            viewer.insert("1." + str(adjusted_start_indx), positions[x])
            viewer.insert("2.8", pos_marker + mini_markers)
        elif x != len(positions):
            position_indx = "1." + str((positions[x] - start_position) + adjusted_start_indx)  # Insert the position
            viewer.insert(position_indx, positions[x])
            position_indx = "2." + str((positions[x] - start_position) + 8) # Insert the marker
            if x != len(positions)-1:
                viewer.insert(position_indx, pos_marker + mini_markers) # Add mini-markers in between positions
            else:
                viewer.insert(position_indx, pos_marker)
    viewer.delete("1." + str(length+20), '1.99999999999') # Remove extra spaces and leave room for the end position
    viewer.delete("2." + str(length+20), '2.99999999999')


# List the fusion genes and add additional info
def add_gene_info(Gene1Name, Gene2Name, Gene1ID, Gene2ID):
    global Gene1Label
    global Gene2Label

    Gene1Link = "https://uswest.ensembl.org/Homo_sapiens/Gene/Summary?g=" + str(Gene1ID)
    Gene2Link = "https://uswest.ensembl.org/Homo_sapiens/Gene/Summary?g=" + str(Gene2ID)

    Gene1Txt.config(state=NORMAL)
    Gene1Txt.insert(INSERT, Gene1Name)
    Gene1Txt.tag_add("here", "1.0", "end")
    Gene1Txt.tag_config("here", background="pale green", foreground="black")
    Gene1Txt.tag_configure("center", justify="center")
    Gene1Txt.tag_add("center", "1.0", "end")
    Gene1Txt.insert(INSERT, "\n\n")
    Gene1Txt.config(state=DISABLED)

    var1 = StringVar()
    Gene1Label = Label(frame3, width=24, height=1, textvariable=var1, bg="pale green", fg="RoyalBlue3", font="Helvetica 11 bold")
    Gene1Label.grid(row=1, column=0, padx=(0,50), pady=(150,0))
    var1.set("More information")
    Gene1Label.bind("<Button-1>", lambda e: show_more_info(Gene1Name, Gene1Link))

    Gene2Txt.config(state=NORMAL)
    Gene2Txt.insert(INSERT, Gene2Name)
    Gene2Txt.tag_add("here", "1.0", "end")
    Gene2Txt.tag_config("here", background="light salmon", foreground="black")
    Gene2Txt.tag_configure("center", justify="center")
    Gene2Txt.tag_add("center", "1.0", "end")
    Gene2Txt.insert(INSERT, "\n\n")
    Gene2Txt.config(state=DISABLED)

    var2 = StringVar()
    Gene2Label = Label(frame3, width=24, height=1, textvariable=var2, bg="light salmon", fg="RoyalBlue3", font="Helvetica 11 bold")
    Gene2Label.grid(row=1, column=2, padx=(50,0), pady=(150,0))
    var2.set("More information")
    Gene2Label.bind("<Button-1>", lambda e: show_more_info(Gene2Name, Gene2Link))

# If fusion is found, display the summary
def add_fusion_summary(Gene1Name, Gene2Name, Gene1Pos, Gene2Pos, Gene1Length, Gene2Length):
    len_diff = Gene1Length - Gene2Length
    intersection = Gene1Pos + len_diff
    Gene1End = str(Gene1Pos + Gene1Length)
    Gene1Chrom = str(samfile.get_chromosome_by_name(Gene1Name[0]))
    Gene2End = str(Gene2Pos + Gene2Length)
    Gene2Chrom = str(samfile.get_chromosome_by_name(Gene2Name[0]))

    SummaryTxt.config(state=NORMAL)
    SummaryTxt.insert(INSERT, Gene1Name[0] + "-" + Gene2Name[0] + "\n\n")
    SummaryTxt.insert(INSERT, "Gene A Read Coordinates:\n")
    SummaryTxt.insert(INSERT, "Chromosome " + str(Gene1Chrom) + ": " + str(Gene1Pos) + " - " + str(Gene1End) +  "\n\n")
    SummaryTxt.insert(INSERT, "Fusion Intersection: " + str(intersection) + "\n\n")
    SummaryTxt.insert(INSERT, "Gene B Read Coordinates:\n")
    SummaryTxt.insert(INSERT, "Chromosome " + str(Gene2Chrom) + ": " + str(Gene2Pos) + " - " + str(Gene2End) +  "\n\n")
    SummaryTxt.tag_configure("center", justify="center")
    SummaryTxt.tag_add("center", "1.0", "end")
    SummaryTxt.config(state=DISABLED)

# Open a window with extra information about the specific gene
def show_more_info(GeneName, GeneLink):
    global infoBox
    subWin = Tk()
    subWin.title(GeneName)

    frame = Frame(subWin)
    frame.pack()

    LinkLabel = Label(frame, width=45, height=1, text="Click for Online Resources", background="light cyan",
                      foreground="blue", font="Helvetica 11 bold")
    LinkLabel.bind("<Button-1>", lambda e: webbrowser.open_new(GeneLink))
    LinkLabel.pack(side=TOP)

    infoBox = Text(frame,  width=45, height=15)

    GeneStart, GeneEnd = samfile.get_positions_by_name(GeneName[0])

    infoBox.insert(INSERT, "ID: " + str(samfile.get_id_by_name(GeneName[0])) + "\n")
    infoBox.insert(INSERT, "Positions: " + str(GeneStart) + " - " + str(GeneEnd) + "\n")
    infoBox.insert(INSERT, "Chromosome: " + str(samfile.get_chromosome_by_name(GeneName[0])) + "\n")

    # Exon regions are identified bases on transcripts
    infoBox.insert(INSERT, "Transcripts: " + "\n")
    trans_ids = samfile.get_transcript_ids_by_name(GeneName[0])
    for t_id in trans_ids:
        t_start, t_end = samfile.get_transcript_positions_by_id(t_id)
        t_name = samfile.get_transcript_name_by_id(t_id)
        infoBox.insert(INSERT, "  " + str(t_name) + "(" + t_id + ")\n  Positions: " + str(t_start) + " - " + str(t_end) + "\n")
        infoBox.insert(INSERT, "    Exon/Intron Boundaries:\n")
        exon_start, exon_end = samfile.get_exon_regions_by_transcript_id(t_id)
        for i in range(len(exon_start)):
            # If the start position does not equal start position of exon, insert intron
            if int(t_start) != exon_start[0] and i == 0:
                infoBox.insert(INSERT, "      I: " + str(t_start) + " - " + str(exon_start[i]-1) + "\n")
            infoBox.insert(INSERT, "      E: " + str(exon_start[i]) + " - " + str(exon_end[i]) + "\n")
            if i != len(exon_start)-1: # If not on the last exon, look ahead to get intron position
                # Insert intron positions as before the start of the next exon and after the end of the current one
                infoBox.insert(INSERT, "      I: " + str(exon_end[i]+1) + " - " + str(exon_start[i+1]-1) + "\n")

    infoBox.pack(side=BOTTOM)
    infoBox.config(state=DISABLED)

    subWin.update()
    subWin.attributes("-topmost", True)
    subWin.resizable(False,False)
    subWin.mainloop()

# This focuses the scroll/text area on the point of fusion genes meeting
def focus_fusion_area(len_diff):
    fusion_indx = "1." + str(len_diff+8)
    fusion_area.see(fusion_indx)

# Exit GUI function
def _quit():
    win.quit()
    win.destroy()
    exit()

# Initialize the main window
def main():
    # Allow the window variables to be accessible outside of the main function (for theme changes, text population, etc.)
    global win
    global frame1
    global frame2
    global frame3
    global frameChr
    global frmChrGen
    global frameRef
    global fusion_area
    global chrom_area
    global genes_in_area
    global refLabel
    global searchlbl1
    global chrlbl1
    global GeneLbl
    global CoordLbl
    global Gene1Txt
    global Gene2Txt
    global SummaryTxt
    global SummaryHeader
    global chrom_header
    global fusion_header
    global Gene1Header
    global Gene2Header
    global scrollbar1
    global sam_open
    global bind_chrom_start
    global bind_chrom_move
    global widget_list # USed for theme changes
    global DEFAULT_THEME

    DEFAULT_THEME = "lightsteelblue1" # The color of the application

    win = Tk()
    win.title("BLUE")
    win.geometry("1000x750")  # Size of window
    win.config(bg=DEFAULT_THEME)

    # Create a menu bar
    menu_bar = Menu(win)
    win.config(menu=menu_bar)

    # Create menu and add menu items
    file_menu = Menu(menu_bar, tearoff=0)
    file_menu.add_command(label='Load from BAM file...', command=_open_sam_file)
    file_menu.add_separator()
    file_menu.add_command(label='Exit', command=_quit)
    menu_bar.add_cascade(label='File', menu=file_menu)
    # Reference human genome menu
    hg_menu = Menu(menu_bar, tearoff=0)
    hg_menu.add_command(label='hg38', command=lambda: _open_genome_file("hg38"))
    menu_bar.add_cascade(label='Ref', menu=hg_menu)
    # Theme options
    theme_menu = Menu(menu_bar, tearoff=0)
    theme_menu.add_command(label="Default", command=lambda: change_theme("lightsteelblue1"))
    theme_menu.add_command(label="Lavender", command=lambda: change_theme("lavender"))
    theme_menu.add_command(label="Rose", command=lambda: change_theme("misty rose"))
    theme_menu.add_command(label="Pinkesque", command=lambda: change_theme("thistle2"))
    theme_menu.add_command(label="Antique", command=lambda: change_theme("antique white"))
    theme_menu.add_command(label="Sky", command=lambda: change_theme("lightskyblue1"))
    theme_menu.add_command(label="Aqua", command=lambda: change_theme("PaleTurquoise2"))
    theme_menu.add_command(label="Snow White", command=lambda: change_theme("snow2"))
    theme_menu.add_command(label="Graydation", command=lambda: change_theme("gray81"))
    menu_bar.add_cascade(label="Theme", menu=theme_menu)

    help_menu = Menu(menu_bar, tearoff=0)
    help_menu.add_command(label="How to...", command=open_how_to)
    help_menu.add_separator()
    help_menu.add_command(label="Report to admin...", command=report_error)
    menu_bar.add_cascade(label="Help", menu=help_menu)

    frameRef = Frame(win)
    frameRef.pack(side=TOP)

    refLabel = Label(frameRef, text='No reference in Use', bg=DEFAULT_THEME)
    refLabel.grid(row=0, column=0)

    frame1 = Frame(win, bg=DEFAULT_THEME)
    frame1.pack()

    # Add a label; enter a gene or a locus
    searchlbl1 = Label(frame1, text='Enter a gene or position:', bg=DEFAULT_THEME)
    searchlbl1.grid(column=3, row=0)

    # Adding a Text box Entry widget for enter a gene or a locus
    search = tk.StringVar()
    search_entered = Entry(frame1, width=12, textvariable=search)
    search_entered.grid(column=4, row=0)

    # Add search button for a gene or a locus
    search_button = Button(frame1, text='Search', command=lambda: search_position(search_entered.get()))
    search_button.grid(column=5, row=0)

    # Add a label with Chromosome drop-down
    chrlbl1 = Label(frame1, text='Chromosome:', bg=DEFAULT_THEME)
    chrlbl1.grid(column=1, row=0)

    # Adding a Text box Entry widget for Chromosome drop-down
    chr = tk.StringVar()
    chr_chosen = ttk.Combobox(frame1, width=12, textvariable=chr)
    chr_chosen['values'] = ("-----Select-----", "chr1", "chr2", "chr3", "chr4", "chr5",
                            "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                            "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
                            "chr18", "chr19", "chr20", "chr21", "chr22", "chrMT", "chrX", "chrY")
    chr_chosen.grid(column=2, row=0)
    chr_chosen.current(0)

    chr_chosen.bind("<<ComboboxSelected>>", load_chromosome)

    frmChrGen = Frame(win, bg=DEFAULT_THEME)
    frmChrGen.pack(side=BOTTOM)

    genes_in_area = Text(frmChrGen, width=50, height=3, wrap=None, bg="snow")
    genes_in_area.grid(row=0, column=1, pady=(0,0))
    genes_in_area.bind("<Enter>", on_enter)
    genes_in_area.tag_configure("center", justify="center")
    genes_in_area.config(state=DISABLED)
    geneTxt = "G\nE\nN\nE\nS"
    GeneLbl = Label(frmChrGen, text=geneTxt, font="Helvetica 12 bold", bg=DEFAULT_THEME, fg="gray27")
    GeneLbl.grid(row=0, column=0, pady=(0,0))
    coordText = "C\nO\nO\nR\nD\n"
    CoordLbl = Label(frmChrGen, text=coordText, font="Helvetica 12 bold", bg=DEFAULT_THEME, fg="gray27")
    CoordLbl.grid(row=0, column=2, pady=(20,0))

    # Initialize the chromosome sequence window
    frameChr = Frame(win, bg=DEFAULT_THEME)
    frameChr.pack(side=BOTTOM)
    chrom_header = Label(frameChr, text="Chrom-view", font="Helvetica 18 bold", bg=DEFAULT_THEME, fg="gray27")
    chrom_header.grid(row=0, column=1, pady=(30,0))
    chrom_area = Text(frameChr, width=107, height=3, wrap=NONE, inactiveselectbackground='white', bg="snow")
    chrom_area.grid(row=1, column=1, pady=(5,0))
    # Add a tag to allow positions markers to indicate between bases
    chrom_area.tag_add("Marker", "1.0", "1.1")
    chrom_area.tag_config("Marker", font=("Georgia", "6", "bold"), foreground="white")
    chrom_area.config(state=DISABLED)
    # Chromosome area events to drag sequence
    chrom_area.bindtags(('Text', 'post-class-bindings', '.' 'all'))
    bind_chrom_start = chrom_area.bind_class('post-class-bindings', '<ButtonPress-1>', move_start)
    bind_chrom_move = chrom_area.bind_class('post-class-bindings', '<B1-Motion>', move_move)
    chrom_area.bind_class('post-class-bindings', '<ButtonRelease-1>', on_release)
    chrom_area.bind_class('post-class-bindings', '<Enter>', on_enter)

    # Initialize the Fusion area
    frame2 = Frame(win, bg=DEFAULT_THEME)
    frame2.pack(side=BOTTOM)
    scrollbar1 = Scrollbar(frame2, orient='horizontal')
    scrollbar1.grid(row=4, column=1, sticky=N + S + E + W, pady=(0, 0))
    fusion_header = Label(frame2, text="Fusion-view", font="Helvetica 18 bold", bg=DEFAULT_THEME, fg="gray27")
    fusion_header.grid(row=2, column=1, pady=(15,0))
    fusion_area = Text(frame2, width=107, height=10, wrap=NONE, xscrollcommand=scrollbar1.set, bg="snow")
    scrollbar1.config(command=fusion_area.xview)
    fusion_area.grid(row=3, column=1, pady=(0,0))

    # Makes it so that text cannot be highlight with cursor
    fusion_area.bindtags((str(fusion_area), str(win), "all"))
    # Add a tag to allow positions markers to indicate between bases
    fusion_area.tag_add("Marker", "1.0", "1.1")
    fusion_area.tag_config("Marker", font=("Georgia","6","bold"), foreground="white")
    fusion_area.config(state=DISABLED)

    # Fusion area events to drag sequences
    fusion_area.bind("<ButtonPress-1>", move_start)
    fusion_area.bind("<B1-Motion>", move_move)
    fusion_area.bind("<ButtonRelease-1>", on_release)
    fusion_area.bind("<Enter>", on_enter)

    frame3 = Frame(win, bg=DEFAULT_THEME)
    frame3.pack(side=BOTTOM)

    # Initialize the fusion boxes
    # Fusion boxes (left & right)
    Gene1Txt = Text(frame3, width=25, height=3, bg="snow")
    Gene1Txt.grid(row=1, column=0, padx=(0,50), pady=(125, 0))
    Gene1Txt.bind("<Enter>", on_enter)
    Gene1Txt.config(state=DISABLED)
    Gene2Txt = Text(frame3, width=25, height=3, bg="snow")
    Gene2Txt.grid(row=1, column=2, padx=(50,0), pady=(125, 0))
    Gene2Txt.bind("<Enter>", on_enter)
    Gene2Txt.config(state=DISABLED)

    Gene1Header = Label(frame3, text="Gene List A", fg="gray17", bg=DEFAULT_THEME, font="Helvetica 12 bold")
    Gene1Header.grid(row=1, column=0, padx=(0,50),pady=(25,0))
    Gene2Header = Label(frame3, text="Gene List B", fg="gray17", bg=DEFAULT_THEME, font="Helvetica 12 bold")
    Gene2Header.grid(row=1, column=2, padx=(50,0),pady=(25,0))

    # Initialize the summary view
    SummaryHeader = Label(frame3, text="Summary", font="Helvetica 12 bold", bg=DEFAULT_THEME, fg="gray27")
    SummaryHeader.grid(row=0, column=1, pady=(0,0))
    SummaryTxt = Text(frame3, width=42, height=10, bg=DEFAULT_THEME, borderwidth=2, relief="groove")
    SummaryTxt.grid(row=1, column=1, pady=(0, 0))
    SummaryTxt.bind("<Enter>", on_enter)
    SummaryTxt.config(state=DISABLED)

    widget_list = [win,frame1,refLabel,searchlbl1,chrlbl1,frmChrGen,GeneLbl,CoordLbl,frameChr,chrom_header,
                   frame2,fusion_header,frame3,Gene1Header,Gene2Header,SummaryHeader,SummaryTxt]

    update_gene_area() # This looks in the chrom-view and updates any visible genes

    win.mainloop()

if __name__=="GUIDesign":
    main()