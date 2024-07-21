# MPI_Recv/MPI_Send versus MPI_Bcast

Implementação de envio de mensagens com MPI_Recv e MPI_Send são menos eficientes: um único processo envia para os outros. Logo, a complexidade é O(n)

Usando o MPI_Bcast temos uma implementação mais eficiente pois o envio de mensagem é feito em árvore. Logo a complexidade é O(log n)

# Mudança do datatype

Nas funções MPI_Gather e MPI_Scatter, temos que informar o datatype do destino e do recebedor.

Um cenário em que isso pode ser útil é quando um processo envia uma matriz indexada por linhas, mas o processo destino faz cálculos indexados por colunas nessa matriz. 

Então, para evitar pulos na memória, enviamos a matriz com um MPI datatype indexada por linhas, e outro processo recebe como um MPI datatype indexada por colunas

Como não existe esses datatypes nativos no MPI, temos que implementá-los.

# Datatypes

MPI_Type_create_struct cria uma estrutura para ser enviada como dado a um processo. 
- array_of_blocklenghts define o tamanho do dado (1 para escalar, por exemplo)
- array_of_displacements define os deslocamentos entre os dados
- array_of_types define os tipos de dados de cada dado

Exemplo. Na regra do trapézio podemos criar um datatype com {a,b,n}. Todos os três são escalares, então iniciamos como `array_of_blocklengths = {1,1,1}`

Ao invés de fazer Bcast com cada um dos três dados separadamente, enviando o tipo de dado customizado input_mpi_t, o compilador só precisa saber o endereço do primeiro elemento na estrutura (que é a_p), e acessar os outros dois utilizando o array de deslocamentos

```c++
MPI_Bcast(a_p, 1, input_mpi_t, 0, MPI_COMM_WORLD)
```
