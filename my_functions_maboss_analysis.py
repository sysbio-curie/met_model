#Libraries
import maboss
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Classification of the phenotypes:
def classify_ss(steady_states, E, M):
    phenotypes={}
    for i in range(len(steady_states)):
        ss=steady_states.iloc[i]
        if(ss[E].any() and ss[M].any()):
            phenotype='P'
        elif (ss[E].any() and not ss[M].any()):
            phenotype='E'
        elif (not ss[E].any() and ss[M].any()):
            phenotype='M'
        else:
            phenotype='U'
        phenotypes[i]= phenotype
    steady_states['Phenotype'] = phenotypes
    return steady_states

def format_phenotypes(steady_states):
    steady_states2=steady_states.copy()

    phenotypes={}

    for i in range(len(steady_states2)):
        if(steady_states2.loc[steady_states2.index[i], 'E'] == 1):
            phenotype='E'
        elif(steady_states2.loc[steady_states2.index[i], 'M'] == 1):
            phenotype='M'
        elif(steady_states2.loc[steady_states2.index[i], 'pEM'] == 1):
            phenotype='pEM'
        else:
            phenotype='U'
        phenotypes[i]= phenotype

    steady_states2['Phenotype'] = phenotypes.values()
    return steady_states2

def get_phenotype(model, nodes_on, output_nodes, nodes_default="off", nodes_random=0, nodes_off=0, summary=0, few_markers_ouput=0):
    on = [0, 1]
    off = [1, 0]
    random = [0.5, 0.5]

    model=model.copy()

    print('Initial conditions \n')
    print('nodes on: ',nodes_on)
    print('rest of nodes: ',nodes_default)

    #set initial conditions
    if (nodes_default == "off"):
        maboss.set_nodes_istate(model, model.network.keys(), off)
    elif(nodes_default == "random"):
        maboss.set_nodes_istate(model, model.network.keys(), random)

    maboss.set_nodes_istate(model, nodes_on, on)

    if (nodes_random):
        maboss.set_nodes_istate(model, nodes_random, random)
        print('nodes random:',nodes_random,'\n')

    if (nodes_off):
        maboss.set_nodes_istate(model, nodes_off, off)
        print('nodes off:',nodes_off,'\n')

    maboss.set_output(model, output_nodes)

    #run simulation
    res = model.run()

    if(summary):

        #get resulting phenotype
        ss=res.get_fptable()
        ss=format_phenotypes(ss)

        if(few_markers_ouput):
            print(ss[["Proba", "Phenotype"]+few_markers_ouput], '\n')
        else:
            print(ss[["Proba","Phenotype"]], '\n')


            print('----------------------------------------------------------')
    return res


def get_phenotype_mutation(model_t, mut_node, mut_type, nodes_on, output_nodes, nodes_default="off", nodes_random=0, nodes_off=0, summary=0, few_markers_ouput=0):
    on = [0, 1]
    off = [1, 0]
    random = [0.5, 0.5]

    model=model_t.copy()

    #print('Initial conditions \n')
    #print('nodes on: ',nodes_on)
    #print('rest of nodes: ',nodes_default)


    #set initial conditions
    if (nodes_default == "off"):
        maboss.set_nodes_istate(model, model.network.keys(), off)
    elif(nodes_default == "random"):
        maboss.set_nodes_istate(model, model.network.keys(), random)

    maboss.set_nodes_istate(model, nodes_on, on)

    if (nodes_random):
        maboss.set_nodes_istate(model, nodes_random, random)
        #print('nodes random:',nodes_random,'\n')

    if (nodes_off):
        maboss.set_nodes_istate(model, nodes_off, off)
        #print('nodes off:',nodes_off,'\n')

    #set output
    maboss.set_output(model, output_nodes)

    #mutate with established initial conditions
    mutant = maboss.copy_and_mutate(model, mut_node, mut_type)
    #run simulation
    res = mutant.run()

    if(summary):

        gene_name=''.join(mut_node)
        print('Mutation', gene_name, mut_type)
        #get resulting phenotype
        ss=res.get_fptable()
        ss=format_phenotypes(ss)

        if(few_markers_ouput):
            print(ss[["Proba", "Phenotype"]+few_markers_ouput], '\n')
        else:
            print(ss[["Proba","Phenotype"]], '\n')


            print('----------------------------------------------------------')
    return res

#TGFB autocrine functions
def tgfb_removal(model, initial_state_on, initial_state_off, dur, dur2, output, col, param):
#dur=30
#output=phenotypes
#initial_state=gene_list_epithelial
    on = [0, 1]
    off = [1, 0]
    random = [0.5, 0.5]

    sim=model.copy()
    sim.param.update(param)
    sim.update_parameters(max_time=dur)
    nodes = list(sim.network.keys())
    sim.network.set_output(output)

    #Define initial state
    initial = sim.copy()
    initial.param["max_time"] = dur
    maboss.set_nodes_istate(initial, initial.network.keys(), random)
    maboss.set_nodes_istate(initial, initial_state_on, on)
    maboss.set_nodes_istate(initial, initial_state_off, off)

    #TGFB signaling starting from the initial state
    sim_short = initial.copy()
    sim_short.network.set_istate('TGFB_L', on)
    #sim_short.update_parameters(max_time=dur)
    sim_short_fig = sim_short.copy()
    res_short_fig = sim_short_fig.run()
    #sim_short.network.set_output(list(sim_short.network.keys()))
    sim_short.network.set_output(nodes)
    res_short = sim_short.run()
    res_short.get_nodes_probtraj()["TGFB_R"].iloc[-1]
    sim_short_fig = sim_short.copy()

    new_istates = change_input(to_istates(res_short.get_states_probtraj(), nodes), "TGFB_L", 0, nodes)

    sim_next = sim_short.copy()
    sim_next.update_parameters(max_time=dur2)
    sim_next.network.set_istate(nodes, new_istates)

    sim_next_fig = sim_next.copy()
    sim_next_fig.network.set_output(output)
    res_next_fig = sim_next_fig.run()

    res_next = sim_next.run()

    res_next.get_nodes_probtraj()["TGFB_R"].iloc[-1]

    #empty_table= res_short_fig.get_states_probtraj().iloc[:1]
    #empty_table=0

    table = res_short_fig.get_states_probtraj()
    table_next = res_next_fig.get_states_probtraj()
    table_next.index = np.array([value + dur - 1 for value in table_next.index.values])

    #col = sns.color_palette('colorblind')

    fig = plt.figure(figsize=(3,2), dpi=200)
    fig.subplots(1)
    pd.concat([table, table_next]).plot(ax=fig.axes[0],color=col)
    #table_next.plot(ax=fig.axes[0],color=col)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8)
    plt.xlabel('Time', fontsize=8)
    plt.ylim(0, 1)
    plt.ylabel('State probability', fontsize=8)
    fig.axes[0].axvspan(0, dur, facecolor='0.9')
    #plt.title("No TGFB autocrine signaling")

def to_bits(state, nodes):
    if state == "<nil>":
        return tuple([0]*len(nodes))
    bits = []
    state_list = state.split(" -- ")
    for node in nodes:
        bits.append(1 if node in state_list else 0)
    return tuple(bits)

def to_istates(table, nodes):
    istates = {}
    for index, value in table.iloc[-1, :].items():
        istates.update({to_bits(index, nodes): value})
    return istates

def change_input(istates, name, value, nodes):
    new_istates = {}
    for bits, proba in istates.items():
        new_bits = list(bits)
        new_bits[nodes.index(name)] = value
        new_tuples = tuple(new_bits)
        if new_tuples not in new_istates.keys():
            new_istates[new_tuples] = proba
        else:
            new_istates[new_tuples] += proba
    return new_istates


def my_pie_chart(maboss_simuation, col):
    fig = plt.figure(figsize=(3,3), dpi=500)
    fig.subplots(1)
    df=maboss_simuation.get_last_states_probtraj()
    #df=df[df >0.01]
    #df=df.dropna(axis=1)

    #if (sum(df.iloc[0])<0.99):
    #    df['Others'] = 1-sum(np.array(df.iloc[0]))

    states=np.array(df.iloc[0])
    names=df.columns
    plt.pie(states, colors=col, autopct='%1.1f%%', startangle=90, normalize=False)
    plt.legend(names, loc="center", bbox_to_anchor=(0.5,-0.1))

def multiple_mutations(initial_condition, mutations, phenotypes):
    #define initial condition
    #initial_condition = initial_condition_t.copy()

    phenotypes_vector={}
    mutation={}
    proportion={}
    i=0

    mutations_keys = list(mutations.keys())
    for mutant in mutations_keys:
        #copy to mutate
        i_condition_mutation= initial_condition.copy()

        #select mutation
        inner_dict=mutations[mutant]

        #iterate over mutations
        for index, (key, value) in enumerate(inner_dict.items()):
            gene=key
            mut_type=value
            i_condition_mutation.mutate(gene, mut_type)

        res = i_condition_mutation.run()
        final_stateprobt=res.get_last_states_probtraj()

        for pheno in phenotypes:
            mutation[i]= mutant
            phenotypes_vector[i]= pheno

            if pheno in list(final_stateprobt.columns):
                proportion[i]= final_stateprobt.iloc[0][pheno]
            else:
                proportion[i]= 0
            i=i+1

    #WT condition that serves as a control
    res= initial_condition.run()
    res.plot_piechart()
    final_stateprobt=res.get_last_states_probtraj()
    n_phenotypes=len(final_stateprobt.columns)

    for pheno in phenotypes:
        mutation[i]= "Wild type"
        phenotypes_vector[i]= pheno

        if pheno in list(final_stateprobt.columns):
            proportion[i]= final_stateprobt.iloc[0][pheno]
        else:
            proportion[i]= 0
        i=i+1

    d = {'mutation': mutation, 'phenotypes': phenotypes_vector, 'proportion':proportion}
    df = pd.DataFrame(data=d)
    return df
