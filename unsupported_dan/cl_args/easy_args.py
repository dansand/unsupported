
def easy_args(sysArgs, _dict):
    """
    Command line args are parsed and values go into the dictionary provided, under key names provided.
    Names are provided in the easyDict format, i.e. dictName.foo, rather than dictName['foo']
    This function will perform different updates, depending on the operator provided:

    * assign new value: =
    * in place multiply: *=
    * in place divide: /=
    * in place add: +=
    * in place subtract: -=



    Parameter
    ---------
    sysArgs : list
        the list of command line args provide by sys.argv (through the sys module)
    _dict:
        a dictionary to update values within

    Examples
    ---------

    >>>python test.py dp.arg3+=24. dp.arg2='Cheryl', dp.arg4=True
    ...


    ...The above code would
    :add 24 to dictionary dp, key 'arg3',
    :replace dictionary dp, key 'arg2', with the string Cheryl
    :replace dictionary dp, key 'arg4', with the Boolean True


    """


    #print(sysArgs)
    for farg in sysArgs:
        #Try to weed out some meaningless args.

        #print(farg)

        if ".py" in farg:
            continue
        if "=" not in farg:
            continue

        try:

            #########
            #Split out the dict name, key name, and value
            #########

            (dicitem,val) = farg.split("=") #Split on equals operator
            (dic,arg) = dicitem.split(".")
            if '*=' in farg:
                (dicitem,val) = farg.split("*=") #If in-place multiplication, split on '*='
                (dic,arg) = dicitem.split(".")
            if '/=' in farg:
                (dicitem,val) = farg.split("/=") #If in-place division, split on '/='
                (dic,arg) = dicitem.split(".")

            if '+=' in farg:
                (dicitem,val) = farg.split("+=") #If in-place addition, split on '+='
                (dic,arg) = dicitem.split(".")
            if '-=' in farg:
                (dicitem,val) = farg.split("-=") #If in-place addition, split on '-='
                (dic,arg) = dicitem.split(".")

            #print(dic,arg,val)

            #########
            #Basic type conversion to float, boolean
            #########

            if val == 'True':
                val = True
            elif val == 'False':     #First check if args are boolean
                val = False
            else:
                try:
                    val = float(val) #next try to convert  to a float,
                except ValueError:
                    pass             #otherwise leave as string

            #########
            #Update the given dictionary
            #########

            try:
                if '*=' in farg:
                        _dict[arg] = _dict[arg]*val #multiply parameter by given factor
                elif '/=' in farg:
                        _dict[arg] = _dict[arg]/float(val) #divide parameter by given value
                elif '+=' in farg:
                        _dict[arg] = _dict[arg]+val #add to parameter given value
                elif '-=' in farg:
                        _dict[arg] = _dict[arg]-val #subtract from parameter given value
                else:
                        _dict[arg] = val    #or reassign parameter by given value

                #
                #else:
                #        _dict[arg] = val    #or reassign parameter by given value

            except:
                    pass

        except:
            pass
