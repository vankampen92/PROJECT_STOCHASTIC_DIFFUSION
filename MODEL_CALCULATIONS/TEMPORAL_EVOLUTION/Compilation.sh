make clean
make MODEL=DIFFUSION_AZTECA_4D BT=PRIORITY_QUEU_SUPER_OPTIMIZATION STO=STOCHASTIC_OPTIMIZATION RE=NON_REUSE_RANDOM_NUMBER
mv DIFFUSION_AZTECA_4D AZTECA_4D_3
make clean
make MODEL=DIFFUSION_AZTECA_4D BT=BINARY_TREE_SUPER_OPTIMIZATION STO=STOCHASTIC_OPTIMIZATION RE=REUSE_RANDOM_NUMBER
mv DIFFUSION_AZTECA_4D AZTECA_4D_22
make clean
make MODEL=DIFFUSION_AZTECA_4D BT=BINARY_TREE_SUPER_OPTIMIZATION STO=STOCHASTIC_OPTIMIZATION RE=NON_REUSE_RANDOM_NUMBER
mv DIFFUSION_AZTECA_4D AZTECA_4D_2
make clean
make MODEL=DIFFUSION_AZTECA_4D BT=NON_BINARY_TREE_SUPER_OPTIMIZATION STO=NON_STOCHASTIC_OPTIMIZATION RE=NON_REUSE_RANDOM_NUMBER
mv DIFFUSION_AZTECA_4D AZTECA_4D_0
make clean
make MODEL=DIFFUSION_AZTECA_4D_0 BT=PRIORITY_QUEU_SUPER_OPTIMIZATION STO=STOCHASTIC_OPTIMIZATION RE=NON_REUSE_RANDOM_NUMBER
mv DIFFUSION_AZTECA_4D_0 AZTECA_4D_03
make clean
make MODEL=DIFFUSION_AZTECA_4D_0 BT=BINARY_TREE_SUPER_OPTIMIZATION STO=STOCHASTIC_OPTIMIZATION RE=REUSE_RANDOM_NUMBER
mv DIFFUSION_AZTECA_4D_0 AZTECA_4D_022
make clean
make MODEL=DIFFUSION_AZTECA_4D_0 BT=BINARY_TREE_SUPER_OPTIMIZATION STO=STOCHASTIC_OPTIMIZATION RE=NON_REUSE_RANDOM_NUMBER
mv DIFFUSION_AZTECA_4D_0 AZTECA_4D_02
make clean
make MODEL=DIFFUSION_AZTECA_4D_0 BT=NON_BINARY_TREE_SUPER_OPTIMIZATION STO=NON_STOCHASTIC_OPTIMIZATION RE=NON_REUSE_RANDOM_NUMBER
mv DIFFUSION_AZTECA_4D_0 AZTECA_4D_00
make clean
