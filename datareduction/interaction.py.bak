alwaysuse = False

################################################################################
# Checks if user wants to do something or not
def yesno(prompt):
    if alwaysuse:
        return True
    ans = ""
    while ans == "":
        ans = raw_input(prompt)
    return (ans[0]=='y' or ans[0]=='Y')

################################################################################
# Sets the default response to yesno
def setalwaysuse(action):
    global alwaysuse
    alwaysuse = action
