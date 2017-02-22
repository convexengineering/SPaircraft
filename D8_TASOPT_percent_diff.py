"""
Method to compute and print percent differences between SP D8 model and TASOPT
"""

def percent_diff(sol):
    """
    Method to compute and print percent differences between SP D8 model and TASOPT

    INPUTS
    ------
    sol: solution from D8 SP model
    """

    #weights to compare
    print "WEIGHT DIFFERENCES"
    print "\n"
    print "Total Fuel Weight % Diff: %s" % compute_diff(mag(sol('W_{f_{total}}')).to('lbf'), 
    print "\n"
    print "Total Aircraft Weight % Diff: %s" % compute_diff(mag(sol('W_{total}')).to('lbf'), 
    print "\n"
    print "Engine Weight % Diff: %s" % compute_diff(mag(sol('W_{engine}')).to('lbf'),
    print "\n"
    print "Tail Weight % Diff: %s" % compute_diff(mag(sol('W_{tail}')).to('lbf'),
    print "\n"
    print "Wing Weight % Diff: %s" % compute_diff(mag(sol('W_{wing}')).to('lbf'),
    print "\n"
    print "Fuselage Weight % Diff: %s" % compute_diff(mag(sol('W_{fuse}')).to('lbf'),
    print "\n"
    print "Payload Weight % Diff: %s" % compute_diff(mag(sol('W_{payload}')).to('lbf'),
    print "\n"
    print "VT Weight % Diff: %s" % compute_diff(mag(sol('W_{VT}')).to('lbf'),
    print "\n"
    print "HT Weight % Diff: %s" % compute_diff(mag(sol('W_{HT}')).to('lbf'),
    print "\n"
    print "Wing Weight % Diff: %s" % compute_diff(mag(sol('W_{wing}')).to('lbf'),

    #wing value to compare
    print "\n\n\n"
    print "WING DIFFERENCES"
    print "\n"
    print "Wing Span % Diff: %s" % compute_diff(mag(sol('b')
    print "\n"
    print "Wing Aspect Ratio % Diff: %s" % compute_diff(mag(sol('AR')
    print "\n"
    print "Wing Area % Diff: %s" % compute_diff(mag(sol('S')
    

    #HT valuees to compare
    print "\n\n\n"
    print "HORIZONTAL TAIL DIFFERENCES"
    print "\n"
    print "HT Aspect Ratio % Diff: %s" % compute_diff(mag(sol('AR_h')
    

    #VT values to compare
    print "\n\n\n"
    print "VERTICAL TAIL DIFFERENCES"
    print "\n"
    print "VT Aspect Ratio % Diff: %s" % compute_diff(mag(sol('A_{vt}')
    print "\n"
    print "TVT Span % Diff: %s" % compute_diff(mag(sol('b_{vt}')
    

def compute_diff(sp, tasopt):
    """
    Method to actually compute the percent difference
    """
    diff = 100*(sp-tasopt)/sp
    
    return diff
