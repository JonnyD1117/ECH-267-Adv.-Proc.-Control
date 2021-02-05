import os 

from pathlib import Path 


if __name__ == '__main__':


    hw_number = int(input("Please Input the HW assignment Number (aka HW#1 or HW#2) NOTE: Just the number!\n"))
    num_probs = int(input("How many problems are required for this HW? \n"))

    root_dir_quest = input("Is current folder of this Python file in the 'Root Directory' of the HW assignment to be made? (Please answer 'y' or 'n'\n)")


    if root_dir_quest == 'y' or root_dir_quest =='Y': 

            root_dir = os.getcwd()
    else: 
            root_dir = input("Please input the PATH to the root directory using (Forward slash '/' convention). Do not include the final /! ")

 
    hw_dir = root_dir + '/HW_'+ str(hw_number)
    latex_build_dir = hw_dir + '/latex_build'
    
    
    img_dir = latex_build_dir + '/images'
    prob_directory = latex_build_dir + '/Problems'
    
    os.mkdir(hw_dir)
    os.mkdir(latex_build_dir)
    os.mkdir(img_dir)

    doc_path = latex_build_dir + "/HW_" + str(hw_number) + "_Solutions_JonathanDorsey.tex"

    Path(doc_path).touch()   


    try: 
        
        os.mkdir(prob_directory)

    except OSError:
        print("Root Problem Folder Creation FAILED")

    for i in range(num_probs): 

        sub_prob_path = prob_directory + '/Prob' + str(i) 
        os.mkdir(sub_prob_path)

        prob_doc = sub_prob_path + '/prob' + str(i) + ".tex"
        Path(prob_doc).touch()

        with open(prob_doc,'a') as prob: 
            prob.write("\\section*{Problem " + str(i) + "}")


        if i == 0:
            with open(doc_path, 'a') as file: 

                file.write("\documentclass[12px]{article} \n")
                file.write("\n")
                file.write("\\usepackage[fleqn]{amsmath} \n")
                file.write("\\usepackage{amssymb} \n")
                file.write("\\usepackage{caption} \n")
                file.write("\\usepackage{subcaption} \n")
                file.write("\\usepackage{graphicx} \n")
                file.write("\\usepackage{cancel} \n")
                file.write("\\usepackage{hyperref} \n")
                file.write("\n")
                file.write("\n")
                file.write("\\hypersetup{colorlinks=true, linkcolor=blue, filecolor=magenta, urlcolor=blue,} \n")
                file.write("\\urlstyle{same} \n")
                file.write("\\graphicspath{ {./images/} } \n")
                file.write("\n")
                file.write("\\newcommand{\R}{\mathbb{R}} \n")
                file.write("\n")
                file.write("\\begin{document} \n")
                file.write("\t \\title{ECH 267 Nonlinear Control Theory \\ Homework \#" + str(hw_number)+ "}\n")
                file.write("\t \\author{Jonathan Dorsey: Department of Mechanical \& Aerospace Engineering} \n")
                file.write("\n")
                file.write("\\maketitle \n")
                file.write("\n")
                file.write("\\begin{center} \n")
                file.write("\t \\section*{Github Repo Hosted at: } \n")
                file.write("\t \\url{https://github.com/JonnyD1117/ECH-267-Adv.-Proc.-Control} \n")
                file.write("\\end{center}\n")
                file.write("\n")

                for j in range(num_probs):
                    file.write("\\include{Problems/Prob" + str(j) + "/prob" + str(j) +"}\n")
                
                file.write("\n")
                file.write("\\end{document} ")

        print("\include{Problem_Files/Prob" + str(i) + "/prob" + str(i)+ "}")







