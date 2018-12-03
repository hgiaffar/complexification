complexification

Our objective is to model how the possibility of un- bounded structural complexification helps any adaptive system to adapt to variable environment. In particular, we are inspired by how i) evolution by natural selection selects genotype-phenotype maps that allow for fast adaptation to changing environments by represent- ing environmental statistics in its structure, and ii) by the ability of cognitive systems to learn to learn, that is, to acquire a general inductive bias that helps switching between tasks of different type and also of different com- plexity. Furthermore, we are interested in how optimal complexification dynamics can be achieved as competing models (i.e. agents) fit the (more and more abundant) environmental data. We therefore impose a simple evo- lutionary dynamics on a set of models, and investigate how different operators of structural variation give rise to different complexification dynamics, while we always make parameter values (as opposed to structural prop- erties) follow random-mutation-and-selection dynamics. In particular, we compare network structures in which i) there is no room for structural complexification, modeled by a network of observed variables and evolving weights between them, ii) structures with fixed complexity, mod- eled by a network of fixed number (and topology?) of hidden variables, and iii) structures with evolvable (hid- den) topology.



functions 

organism: class containing constructor function, variational operators and fitness functions

RouletteWheel: fitness proportional selection

evo_dynamics_step: forward step applying - k-cycle Moran model


