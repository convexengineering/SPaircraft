def geomFlag(aircraft):
    global D80, D82, D82, D82_73eng, D8_eng_wing, D8big, b737800, b777300ER, optimal737, \
        optimalD8, Mo8D8, M08_D8_eng_wing, M072_737, D8fam, D8_no_BLI, \
        M08D8_noBLI, optimal777, D8big_eng_wing, multimission, \
        D8bigfam, optimalRJ, RJfam, smallD8, smallD8_no_BLI, smallD8_eng_wing, D12
    global wingengine, rearengine, doublebubble, tube, piHT, conventional
    if aircraft == 'D80':
        D80 = True;
        rearengine = True;
        BLI = True;
        piHT = True;
        doublebubble = True;
        eng = 3;
    if aircraft == 'D82':
        D82 = True;
        rearengine = True;
        BLI = True;
        piHT = True;
        doublebubble = True;
        eng = 3;
    if aircraft == 'D82_73eng':
        D82_73eng = True;
        rearengine = True;
        BLI = True;
        piHT = True;
        doublebubble = True;
        eng = 1;
    if aircraft == 'D8_eng_wing':
        D8_eng_wing = True;
        wingengine = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'D8big':
        D8big = True;
        rearengine = True;
        BLI = True;
        piHT = True;
        doublebubble = True;
        eng = 4;
    if aircraft == 'D8big_no_BLI':
        D8big_no_BLI = True;
        rearengine = True;
        piHT = True;
        doublebubble = True;
        eng = 4;
    if aircraft == 'D8big_eng_wing':
        D8big_eng_wing = True;
        wingengine = True;
        piHT = True;
        doublebubble = True;
        eng = 4;
    if aircraft == 'D8big_M072':
        D8big = True
        D8big_M072 = True;
        rearengine = True;
        BLI = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'D8big_no_BLI_M072':
        D8big_no_BLI = True
        D8big_no_BLI_M072 = True;
        rearengine = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'D8big_eng_wing_M072':
        D8big_eng_wing = True
        D8big_eng_wing_M072 = True;
        wingengine = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'b737800':
        b737800 = True;
        conventional = True;
        eng = 1;
    if aircraft == 'b777300ER':
        b777300ER = True;
        conventional = True;
        eng = 4;
    if aircraft == 'optimal737':
        optimal737 = True;
        conventional = True
    if aircraft == 'optimalD8':
        optimalD8 = True;
        rearengine = True;
        BLI = True;
        piHT = True;
        doublebubble = True;
        eng = 3;
    if aircraft == 'optimal777':
        optimal777 = True;
        conventional = True;
        eng = 4;
    if aircraft == 'optimal777_M08':
        optimal777 = True
        optimal777_M08 = True;
        conventional = True
    if aircraft == 'optimal777_M072':
        optimal777 = True
        optimal777_M072 = True;
        conventional = True
    if aircraft == 'M08D8':
        M08D8 = True;
        rearengine = True;
        BLI = True;
        piHT = True;
        doublebubble = True;
        eng = 3;
    if aircraft == 'M08D8_noBLI':
        M08D8_noBLI = True;
        rearengine = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'M08_D8_eng_wing':
        M08_D8_eng_wing = True;
        wingengine = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'M072_737':
        M072_737 = True;
        conventional = True
    if aircraft == 'D8_no_BLI':
        D8_no_BLI = True;
        rearengine = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'D8big_M072':
        D8big = True
        D8big_M072 = True;
        rearengine = True;
        BLI = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'D8big_M08':
        D8big = True
        D8big_M08 = True;
        rearengine = True;
        BLI = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'optimalRJ':
        optimalRJ = True;
        conventional = True
    if aircraft == 'smallD8':
        smallD8 = True;
        rearengine = True;
        BLI = True;
        piHT = True;
        doublebubble = True;
        eng = 3;
    if aircraft == 'smallD8_eng_wing':
        smallD8_eng_wing = True;
        wingengine = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'smallD8_no_BLI':
        smallD8_no_BLI = True;
        rearengine = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'smallD8_M08_no_BLI':
        smallD8_no_BLI = True
        smallD8_M08_no_BLI = True;
        rearengine = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'smallD8_M08':
        smallD8 = True
        smallD8_M08 = True;
        rearengine = True;
        BLI = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'smallD8_M08_eng_wing':
        smallD8_eng_wing = True
        smallD8_M08_eng_wing = True;
        wingengine = True;
        piHT = True;
        doublebubble = True;
    if aircraft == 'D12':
        D12 = True;
        rearengine = True;
        BLI = True;
        piHT = True;
        doublebubble = True;
        eng = 4;
