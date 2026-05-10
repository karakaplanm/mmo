
import os
import glob


def citation():
    print("********************************************************************************")
    print("please cite this plugin in your work The doi is: 10.5281/zenodo.6821116")
    print("get other citation information in my github page https://github.com/karthikn130")
    print("********************************************************************************")
    
def help():
    print("***************************************************************************")
    print("Version 2.0")
    print("This tool shows the docking scores from ")
    print("                1. Autodock dlg file located in subfolders")
    print("                2. Autodock vina output log file located in subfolders")
    print("                3. Autodock vina output pdbqt file located in subfolders")
    print("To view all docking scores from working directory or Main docking folder")
    print("Enter the Main docking folder location and press submit")
    print("Dock folder - working directory of autodock")
    print("Main docking folder - Folder which contain all docking folders")
    print("***************************************************************************")
    print("Thank you for using this plugin")
    print("For feedback write to karthikn130@gmail.com")
    print("For more detailed help visit PymolBiomolecules Youtube Chennel")


def check_list(list):
    #check if a list is not empty
    # if list is not empty return true
    # if list is empty return false
    if list:
        return True
    else:
        return False


def score_sort(scores):
    # sorts the scores and returns the best score

    # check if the list is empty
    if check_list(scores) == False:
        return "No scores"
    scores.sort()
    best_conf = scores[0]
    return best_conf



def check_extension(main_folder, extension):
    #creates a list of files with that ext in folder and returns the list

    #check if main_folder is a file
    if os.path.isfile(main_folder):
        if main_folder.endswith(extension):
            return [main_folder]

    # get files in main_folder
    folder_ext = os.path.join(main_folder, extension)
    file_list = []
    file_list = glob.glob(folder_ext)
    

    # get files in subfolders also
    sub_fol = os.listdir(main_folder)
    for fol in sub_fol:
        if os.path.isdir(os.path.join(main_folder, fol)) == True:
            folder_ext = os.path.join(main_folder, fol, extension)
            sub_list = glob.glob(folder_ext)
            if check_list(sub_list) == True:
                file_list.extend(sub_list)
  
    return file_list



def parse_dlg(main_folder):
    # this function parses the dlg files in folder and its subfolders
    dlg_list = check_extension(main_folder, "*.dlg")
    for dlg in dlg_list:
        search = "DOCKED: USER    Estimated Free Energy of Binding"
        fhandle = open(dlg, encoding='utf-8')
        lines = fhandle.readlines()
        scores = []
        for line in lines:
            if line.startswith(search):
                word = line.split()
                score = word[8]
                #list = score.split(sep = "-")
                #number = list[1]
                number = float(score)
                scores.append(number)
        best_conf = score_sort(scores)
        print(dlg, best_conf)
    citation()



def parse_vina_log(vina_log_folder):
    # this function parses the vina output log file
    log_list = check_extension(vina_log_folder, "*.txt")
    for log in log_list:
        # print name of the text file 
        print(log)
        search = "   1"
        ligand_name = "Output will be"
        vina_log_file = open(log, "r")
        lines = vina_log_file.readlines()
        for line in lines:
            if line.startswith(ligand_name):
                word = line.split()
                ligand = word[3]

            if line.startswith(search):
                word = line.split()
                score = word[1]
                print(ligand, score)
    citation()


def parse_vina_pdbqt(vina_pdbqt_fol):
    # takes all vina output pdbqt files in a folder and prints the best score
    pdbqt_list = check_extension(vina_pdbqt_fol, "*.pdbqt")
    for pdbqt in pdbqt_list:
        search = "REMARK VINA RESULT:"
        vina_pdbqt_file = open(pdbqt, "r")
        lines = vina_pdbqt_file.readlines()
        scores = []
        for line in lines:
            if line.startswith(search):
                word = line.split()
                score = word[3]
                scores.append(float(score))
        best_conf = score_sort(scores)
        # print only if best_conf is not "No scores"
        if best_conf != "No scores":
            print(pdbqt, best_conf)
    citation()



# while testing the plugin uncomment the following lines
# run if file is main
# this file can be run independently after this point
# if __name__ == "__main__":
#     print("This is a module")
#     print("This plugin is intended to be used in pymol")
#     citation()