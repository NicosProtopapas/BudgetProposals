from gurobipy import *
import sys, getopt

def usage():
    print( 'loss_calculation.py --help --optimal --zero --nonzero')
    print()
    print("-h, --help: shows this help message" )
    print("-o, --optimal: only programs with optimal solutions printed")
    print("-n, --nonzero: computes the upper bound only for the case x_1>0,x_2>0,x_3>0")
    print("-z, --zero: computes the upper bound only for the case x_3=0")

def pattern_to_string(pattern):
    """
        A function to convert sign patterns to strings.
    """
    s=""
    for i in pattern:
        if i==1:
            s=s+'+'
        if i==-1:
            s=s+'-'
    return s;

def print_helper(results,suppress_infeasible=False):

    print ("{:<10} {:<10} {:<10} {:<15} {:<30} {:<30} {:<30}".format('t','Sign Pattern','Phantom Pattern','Status','Best Feasible Sol.','Best Upper Bound','Difference'))
    print()
    counter=1;
    for k,v in results.items():
        T,sign,pattern,status,loss,upper,diff=v
        if(suppress_infeasible==True and status=='INFEASIBLE'):
            counter=counter+1;
            continue
        if(status=='INFEASIBLE'):
            print ("{:<10} {:<10} {:<10} {:<15} {:<30} {:<30} {:<30}".format(T,sign,pattern,status,loss,upper,'-'))
        else:
            print ("{:<10} {:<10} {:<10} {:<15} {:<30} {:<30} {:<30}".format(T,sign,pattern,status,loss,upper,diff))
    if(suppress_infeasible):
        print(counter,"infeasible programs are hidden")

def main():
    suppress_infeasible=False
    only_zero=False
    only_nonzero=False
    try:
        opts, args = getopt.getopt(sys.argv[1:], "honz", ["help", "optimal", "nonzero","zero"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit
    print(opts)
    for o,a in opts:
        if o in ("-h","--help"):
            usage()
            sys.exit()
        elif o in ("-o","--optimal"):
            suppress_infeasible=True
        elif o in ("-n","--nonzero"):
            only_nonzero=True
        elif o in ("-z","--zero"):
            only_zero=True
        else:
            usage()
            assert False, "unhandled option"

    #solver(suppress_infeasible,print_location)

    if(only_nonzero):
        results=solver()
        print("Printing the results for x_1>0,x_2>0,x_3>0")
        print_helper(results,suppress_infeasible)
        print()
    elif(only_zero):
        results_special=special_case();
        print("Printing the results for x_1>0,x_2>0,x_3=0")
        print_helper(results_special,suppress_infeasible)
        print()
    else:
        results_special=special_case();
        results=solver()
        print("Printing the results for x_1>0,x_2>0,x_3>0")
        print_helper(results,suppress_infeasible)
        print()
        print("Printing the results for x_1>0,x_2>0,x_3=0")
        print_helper(results_special,suppress_infeasible)
        print()



def solver():

    """
        A function that creates and runs a series of Quadratic Programs.
        These programs finds the three-type profile with the maximum loss.
        We are in the case where  x_1>0,x_2>0,x_3>0.

        Each Quadratic Program is build based on 3 levels.
        1) If t >=1/2 ro t<=1/2
        2) 2 sign patterns -- the signs of \bar{v}_j-x_j
        3) phantom patterns: Whether x_j is between two red phantoms or two black
        (It cannot be between a black and a red, or vice versa for three-type profiles)

    """

    sign_patterns={1:[1,-1,-1],2:[1,1,-1]}
    side_of_t={1:"t>=1/2",2:"t<=1/2"}
    results=dict()

    no_ph_patterns=3**3; # number of phantom patterns
    temp_phantom_patterns=[x for x in itertools.product(['(b,b)','(b,r)','(r,r)'],repeat=3)]
    phantom_patterns=dict(zip(range(1,no_ph_patterns+1),temp_phantom_patterns))

    projects=[1,2,3]
    distinct_pairs=[(i,j) for (i,j) in itertools.product(projects,projects) if i!=j]

    model=dict();

    # handlers for the variables;
    t=dict();
    a=dict();
    b=dict();
    x=dict();
    A=dict();
    B=dict();

    print("Initializing 32 Quadradic Programs with Quadradic Constraints")

    """
        We are using the tuple (Wp,Wt,i) to enumerate our programs according to
        the sign pattern, if t>=1/2 or not, and the phantom pattern.
    """

    for Wp in range(1,2+1):
        for Wt in range(1,2+1):
            for i in range(1,no_ph_patterns+1):
                model_name=pattern_to_string(sign_patterns[Wp])+':'+''.join(phantom_patterns[i])+':'+side_of_t[Wt]
                working_tuple=(Wp,Wt,i)

                model[working_tuple]=Model(model_name) # create a new model
                model[working_tuple].reset() # clear any previous instances
                model[working_tuple].params.NonConvex = 2 # All of our programs are non-convex.
                model[working_tuple].params.MIPGap=1e-5 # We set the preffered error tolerence

                if Wt==1:
                    t[working_tuple] = model[working_tuple].addVar(vtype=GRB.CONTINUOUS,name='t',lb=1/2,ub=1)
                if Wt==2:
                    t[working_tuple] = model[working_tuple].addVar(vtype=GRB.CONTINUOUS,name='t',lb=0,ub=1/2)
                a[working_tuple] = model[working_tuple].addVars(projects,vtype=GRB.CONTINUOUS,name='a',lb=0,ub=1)
                x[working_tuple] = model[working_tuple].addVars(projects,vtype=GRB.CONTINUOUS,name='x',lb=0,ub=1)
                b[working_tuple] = model[working_tuple].addVars(distinct_pairs,vtype=GRB.CONTINUOUS,name='b',lb=0,ub=1)
                model[working_tuple].update()

    """
        Define Objectives: For the sign pattern (+--) and t >=1, i.e. when Wp==1
        and Wt==1 we maximize th get the maximum loss. For the rest Programs,
        we check if any one yields loss greater than 2/3, bu adding the objective
        function as a contraint.
    """


    for Wp in range(1,2+1):
        for Wt in range(1,2+1):
            for i in range(1,no_ph_patterns+1):
                working_tuple=(Wp,Wt,i)
                if (Wp==1 and Wt==1):

                    A[working_tuple]=quicksum(a[working_tuple][j] for j in projects)
                    B[working_tuple]=quicksum(b[working_tuple][j] for j in distinct_pairs )
                    v=dict();
                    for j in projects:
                        v[j]=a[working_tuple][j]+(1-A[working_tuple]-B[working_tuple])*x[working_tuple][j] + quicksum((1-x[working_tuple][k])*b[working_tuple][k,j] for k in projects if k!=j) + x[working_tuple][j]*quicksum(b[working_tuple][j,k] for k in projects if k!=j)

                    signs=sign_patterns[Wp]
                    loss=quicksum(signs[j-1]*(v[j]-x[working_tuple][j]) for j in projects)

                    model[working_tuple].setObjective(loss,GRB.MAXIMIZE)
                    model[working_tuple].update()
                else:
                    A[working_tuple]=quicksum(a[working_tuple][j] for j in projects)
                    B[working_tuple]=quicksum(b[working_tuple][j] for j in distinct_pairs )
                    v=dict();
                    for j in projects:
                        v[j]=a[working_tuple][j]+(1-A[working_tuple]-B[working_tuple])*x[working_tuple][j] + quicksum((1-x[working_tuple][k])*b[working_tuple][k,j] for k in projects if k!=j) + x[working_tuple][j]*quicksum(b[working_tuple][j,k] for k in projects if k!=j)

                    signs=sign_patterns[Wp]
                    loss=quicksum(signs[j-1]*(v[j]-x[working_tuple][j]) for j in projects)

                    model[working_tuple].setObjective(1,GRB.MAXIMIZE)
                    model[working_tuple].addConstr(loss>=2/3)
                    model[working_tuple].update()

    """
    We add the constraints:
    """

    for Wp in range(1,2+1):
        for Wt in range(1,2+1):
            for i in range(1,no_ph_patterns+1):
                working_tuple=(Wp,Wt,i)
                phantoms=phantom_patterns[i]
                signs=sign_patterns[Wp]

                z=dict()
                q=dict()
                for j in projects:
                    z[j]= a[working_tuple][j] + quicksum(b[working_tuple][k,j] for k in projects if k!=j)
                    q[j]= quicksum(b[working_tuple][j,k] for k in projects if k!=j)

                # constraints added to all models
                # feasible solution and A+B<=1
                model[working_tuple].addConstr(quicksum(x[working_tuple][j] for j in projects)==1)
                model[working_tuple].addConstr(A[working_tuple]+B[working_tuple]<=1)

                if Wt==1: # when t>1/2
                    for j in projects:
                        # ['(b,b)','(b,r)','(r,r)']
                        # recall that C=1-A[working_tuple]-B[working_tuple]
                        if phantoms[j-1]=='(b,b)':
                            # The median is between two black phantoms
                            model[working_tuple].addConstr( x[working_tuple][j]>= z[j]*(2*t[working_tuple]-1))
                            model[working_tuple].addConstr( x[working_tuple][j]<=(1-A[working_tuple]-B[working_tuple]+z[j]+q[j])*(2*t[working_tuple]-1))
                            model[working_tuple].addConstr( z[j]<=1/2)
                            model[working_tuple].addConstr((1-A[working_tuple]-B[working_tuple]+z[j]+q[j])<=1/2)
                        if phantoms[j-1]=='(b,r)':
                            # The median is between a black (lower bound) and a red phantom (upper bound)
                            model[working_tuple].addConstr( x[working_tuple][j]>= z[j]*(2*t[working_tuple]-1))
                            model[working_tuple].addConstr( x[working_tuple][j]<=(1-A[working_tuple]-B[working_tuple]+z[j] + q[j] )*(3-2*t[working_tuple])-2+2*t[working_tuple])
                            model[working_tuple].addConstr( z[j]<=1/2)
                            model[working_tuple].addConstr( (1-A[working_tuple]-B[working_tuple]+z[j] + q[j] )>=1/2)
                        if phantoms[j-1]=='(r,r)':
                            # The median is between two red phantoms
                            model[working_tuple].addConstr( x[working_tuple][j]>=z[j]*(3-2*t[working_tuple])-2+2*t[working_tuple])
                            model[working_tuple].addConstr( x[working_tuple][j]<=(1-A[working_tuple]-B[working_tuple]+z[j] + q[j] )*(3-2*t[working_tuple])-2+2*t[working_tuple])
                            model[working_tuple].addConstr( z[j]>=1/2)
                            model[working_tuple].addConstr( (1-A[working_tuple]-B[working_tuple]+z[j] + q[j] )>=1/2)
                if Wt==2: # when t<=1/2
                    for j in projects:
                        # The median is between two black phantoms
                        if phantoms[j-1]=='(b,b)':
                            model[working_tuple].addConstr(x[working_tuple][j]>=0)
                            model[working_tuple].addConstr(x[working_tuple][j]<=0)
                            model[working_tuple].addConstr(z[j] <=1/2)
                            model[working_tuple].addConstr( (1-A[working_tuple]-B[working_tuple]+z[j]+q[j]) <=1/2 )
                        if phantoms[j-1]=='(b,r)' :
                        # The median is between a black (lower bound) and a red phantom (upper bound)
                            model[working_tuple].addConstr(x[working_tuple][j]>=0)
                            model[working_tuple].addConstr(x[working_tuple][j]<=(1-A[working_tuple]-B[working_tuple]+z[j])*(4*t[working_tuple])-2*t[working_tuple])
                            model[working_tuple].addConstr(z[j]<=1/2)
                            model[working_tuple].addConstr(( 1-A[working_tuple]-B[working_tuple]+z[j]+q[j] )>=1/2)
                        if phantoms[j-1]=='(r,r)' :
                        # The median is between a black (lower bound) and a red phantom (upper bound)
                            model[working_tuple].addConstr(x[working_tuple][j]>=z[j]*(4*t[working_tuple])-2*t[working_tuple])
                            model[working_tuple].addConstr(x[working_tuple][j]<=(1-A[working_tuple]-B[working_tuple]+z[j])*(4*t[working_tuple])-2*t[working_tuple])
                            model[working_tuple].addConstr(z[j]>=1/2)
                            model[working_tuple].addConstr((1-A[working_tuple]-B[working_tuple]+z[j]+q[j])>=1/2)
                model[working_tuple].update()

    print("All models created and intialized. Run the programs")

    for Wp in range(1,2+1):
        for Wt in range(1,2+1):
            for k in range(1,no_ph_patterns+1):
                model[(Wp,Wt,k)].optimize()

    """
    Collect the results in a presentable enviroment
    """
    for Wp in range(1,2+1):
        for Wt in range(1,2+1):
            for k in range(1,no_ph_patterns+1):
                if model[(Wp,Wt,k)].status==2:
                    diff=model[(1,1,k)].objBound-model[(1,1,k)].objVal;
                    results[(Wp,Wt,k)]=(side_of_t[Wt],pattern_to_string(sign_patterns[Wp]),''.join(phantom_patterns[k]),"OPTIMAL",model[(1,1,k)].objVal,model[(1,1,k)].objBound,diff)
                else:
                    results[(Wp,Wt,k)]=(side_of_t[Wt],pattern_to_string(sign_patterns[Wp]),''.join(phantom_patterns[k]), "INFEASIBLE",'-','-','-' )

    return results

def special_case():
    """
        A function that creates and runs a series of Quadratic Programs.
        These programs finds the three-type profile with the maximum loss.
        We are in the special case where  x_1>0,x_2>0,x_3=0.

        Only for $t<1/2$ (Case $t>1/2$ is descarded theoretically )
        Two sign patterns: (-,-,+),(-,+,+)
        phantom patterns are different here: the median is no always between
        two phantoms of the same color. Instead, is :
        case (1): either between two red phantoms or
        case (2): is lower bounded by a black and upper bounded by a red.

        (special_flag,Wp,Wt,phantoms)
    """

    sign_patterns_special={1:[-1,1,1],2:[-1,-1,1]}
    phantom_patterns={1:"rr-rr",2:"rb-rr",3:"rr-rb",4:"rb-rb"}
    side_of_t={1:"t>=1/2",2:"t<=1/2"}

    results=dict()

    projects=[1,2,3]

    distinct_pairs=[(i,j) for (i,j) in itertools.product(projects,projects) if i!=j]

    model=dict();

    # handlers for the variables;
    t=dict();
    a=dict();
    b=dict();
    x=dict();
    A=dict();
    B=dict();

    special=1 # special flag is on;
    Wt=2;
    for Wp in range(1,2+1):
        for i in range(1,4+1):
            model_name=pattern_to_string( sign_patterns_special[Wp] )+':'+''.join(phantom_patterns[i])+':'+side_of_t[Wt]
            working_tuple=(special,Wp,Wt,i)

            model[working_tuple]=Model('model_name') # create a new model
            model[working_tuple].reset() # clear any previous instances
            model[working_tuple].params.NonConvex = 2 # All of our programs are non-convex.
            model[working_tuple].params.MIPGap=1e-5 # We set the preffered error tolerence


            if Wt==1:
                t[working_tuple] = model[working_tuple].addVar(vtype=GRB.CONTINUOUS,name='t',lb=1/2,ub=1)
            if Wt==2:
                t[working_tuple] = model[working_tuple].addVar(vtype=GRB.CONTINUOUS,name='t',lb=0,ub=1/2)
            a[working_tuple] = model[working_tuple].addVars(projects,vtype=GRB.CONTINUOUS,name='a',lb=0,ub=1)
            x[working_tuple] = model[working_tuple].addVars(projects,vtype=GRB.CONTINUOUS,name='x',lb=0,ub=1)
            b[working_tuple] = model[working_tuple].addVars(distinct_pairs,vtype=GRB.CONTINUOUS,name='b',lb=0,ub=1)

    """
    Add the objective as a contraint
    """
    for Wp in range(1,2+1):
        for i in range(1,4+1):
            working_tuple=(special,Wp,Wt,i)
            A[working_tuple]=quicksum(a[working_tuple][j] for j in projects)
            B[working_tuple]=quicksum(b[working_tuple][j] for j in distinct_pairs )

            v=dict();
            for j in projects:
                v[j]=a[working_tuple][j]+(1-A[working_tuple]-B[working_tuple])*x[working_tuple][j] + quicksum((1-x[working_tuple][k])*b[working_tuple][k,j] for k in projects if k!=j) + x[working_tuple][j]*quicksum(b[working_tuple][j,k] for k in projects if k!=j)

            signs=sign_patterns_special[Wp]
            loss=quicksum(signs[j-1]*(v[j]-x[working_tuple][j]) for j in projects)

            model[working_tuple].setObjective(loss,GRB.MAXIMIZE)
            #model[working_tuple].addConstr(loss>=2/3)

    for Wp in range(1,2+1):
        for i in range(1,4+1):
            working_tuple=(special,Wp,Wt,i)

            z=dict()
            model[working_tuple].addConstr( quicksum( x[working_tuple][j] for j in projects) ==1)
            model[working_tuple].addConstr(A[working_tuple]+B[working_tuple]<=1)
            model[working_tuple].addConstr(x[working_tuple][3]==0)

            z[1]=1-a[working_tuple][2]-a[working_tuple][3]-b[working_tuple][2,3]-b[working_tuple][3,2];
            z[2]=1-a[working_tuple][1]-a[working_tuple][3]-b[working_tuple][1,3]-b[working_tuple][3,1];

            if i==1:
                # rr-rr
                # project 1
                model[working_tuple].addConstr(x[working_tuple][1] >= (a[working_tuple][1] + b[working_tuple][1,3])*4*t[working_tuple] -2*t[working_tuple]  )
                model[working_tuple].addConstr(a[working_tuple][1]+ b[working_tuple][1,3]>=1/2)
                model[working_tuple].addConstr(x[working_tuple][1] <= z[1]*4*t[working_tuple] -2*t[working_tuple]  )
                model[working_tuple].addConstr(z[1]>=1/2)
                # project 2
                model[working_tuple].addConstr(x[working_tuple][2] >= (a[working_tuple][2]+ b[working_tuple][3,2])*4*t[working_tuple] -2*t[working_tuple]  )
                model[working_tuple].addConstr(a[working_tuple][2]+ b[working_tuple][3,2]>=1/2)
                model[working_tuple].addConstr(x[working_tuple][2] <= z[2]*4*t[working_tuple] -2*t[working_tuple]  )
                model[working_tuple].addConstr(z[2]>=1/2)
            if i==2:
                # rb-rr
                # project 1
                model[working_tuple].addConstr(x[working_tuple][1] >= 0  )
                model[working_tuple].addConstr(a[working_tuple][1]<=1/2)
                model[working_tuple].addConstr(x[working_tuple][1] <= z[1]*4*t[working_tuple] -2*t[working_tuple]  )
                model[working_tuple].addConstr(z[1]>=1/2)
                # project 2
                model[working_tuple].addConstr(x[working_tuple][2] >= a[working_tuple][2]*4*t[working_tuple] -2*t[working_tuple]  )
                model[working_tuple].addConstr(a[working_tuple][2]>=1/2)
                model[working_tuple].addConstr(x[working_tuple][2] <= z[2]*4*t[working_tuple] -2*t[working_tuple]  )
                model[working_tuple].addConstr(z[2]>=1/2)
            if i==3:
                # rr-rb
                # project 1
                model[working_tuple].addConstr(x[working_tuple][1] >= a[working_tuple][1]*4*t[working_tuple] -2*t[working_tuple]  )
                model[working_tuple].addConstr(a[working_tuple][1]>=1/2)
                model[working_tuple].addConstr(x[working_tuple][1] <= z[1]*4*t[working_tuple] -2*t[working_tuple]  )
                model[working_tuple].addConstr(z[1]>=1/2)
                # project 2
                model[working_tuple].addConstr(x[working_tuple][2] >= 0  )
                model[working_tuple].addConstr(a[working_tuple][2]<=1/2)
                model[working_tuple].addConstr(x[working_tuple][2] <= z[2]*4*t[working_tuple] -2*t[working_tuple]  )
                model[working_tuple].addConstr(z[2]>=1/2)
            if i==4:
                # rb-rb
                # project 1
                model[working_tuple].addConstr(x[working_tuple][1] >= 0  )
                model[working_tuple].addConstr(a[working_tuple][1]<=1/2)
                model[working_tuple].addConstr(x[working_tuple][1] <= z[1]*4*t[working_tuple] -2*t[working_tuple]  )
                model[working_tuple].addConstr(z[1]>=1/2)
                # project 2
                model[working_tuple].addConstr(x[working_tuple][2] >= 0  )
                model[working_tuple].addConstr(a[working_tuple][2]<=1/2)
                model[working_tuple].addConstr(x[working_tuple][2] <= z[2]*4*t[working_tuple] -2*t[working_tuple]  )
                model[working_tuple].addConstr(z[2]>=1/2)


            model[working_tuple].update()

            print("All models created and intialized. Run the programs")

    for Wp in range(1,2+1):
        for k in range(1,4+1):
            model[(special,Wp,Wt,k)].optimize()

    """
    Collect the results for printing
    """
    for Wp in range(1,2+1):
        for k in range(1,4+1):
            if model[(special,Wp,Wt,k)].status==2:
                diff=model[(special,Wp,Wt,k)].objBound-model[(special,Wp,Wt,k)].objVal;
                results[(special,Wp,Wt,k)]=(side_of_t[Wt],pattern_to_string(sign_patterns_special[Wp]),''.join(phantom_patterns[k]),"OPTIMAL",model[(special,Wp,Wt,k)].objVal,model[(special,Wp,Wt,k)].objBound,(diff))
            else:
                results[(special,Wp,Wt,k)]=(side_of_t[Wt],pattern_to_string(sign_patterns_special[Wp]),''.join(phantom_patterns[k]), "INFEASIBLE",'-','-','-' )

    return(results)
    # End of main

if __name__ == "__main__":
    main()
