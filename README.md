##  DINATEMP
**Simulação de moléculas de Argônio utilizando scilab em ambiente paralelo**

O DINATEMP é um programa criado com linguagem Scilab como ambiente programático e sistema operacional *Linux Debian 8.* 

No sistema criado ocorrem interações em escala molecular descritas através de suas  posições e momentos angulares. 

Através de simulações da dinâmica molecular desses sistemas e com auxilio de métodos 
teóricos pode-se obter a descrição das interações intermoleculares e a definição dos 
núcleos atômicos que interagem sob um potencial.

As condições iniciais da simulação serão geradas da seguinte forma: 
* Criação de um modelo inicial com a geometria aproximada do sistema para os elementos envolvidos dentro de uma caixa de simulação;
* Execução os parâmetros iniciais da simulação, deve-se anexar ao sistema um fator de transferência de energia do/para reservatório de modo a se alcançar o estado de equilíbrio o mais rápido possível;
* Estabelecer ao programa configuração inicial do sistema antes da simulação;
* Localização das paredes da caixa de simulação (por eixo cartesiano);
* Definir a temperatura de equilíbrio do sistema.

As simulações são efetuadas e as trajetórias definidas pela integração das equações 
de **Newton via Runge-Kutta**
Há uma variação de energia durante o tempo de simulação. 

Para iniciar a simulação nessa parte utilizou-se os dados do **Cambridge** **Cluster**  **Database**.

```Cambridge Cluster Database é um banco de dados com N-posições e  energias de clusters diversos  com base em experimentos teóricos ``` 

Também nas condições iniciais, temos os momentos angulares *(px, py , pz)* 
Os vetores de entrada compostos pelas 6 coordenadas atômicas de posição 
e momentos *(x, y, z, px, py , pz)*.
O tempo de análise para as simulações são em pico-segundos, o que facilita a resolução 
das equações ordinárias.

``` Passos de integração menores levam a uma melhor integração, porém demandam tempo computacional maior ``` 

O cálculo da simulação em tempo real se dá através do produto entre o passo e a temperatura  final da simulação (TS), em unidades reduzidas. 
É importantíssimo que haja uma conversão do tempo de simulação para o tempo real, dando origem a um vetor de tempo real (TR) descrito em forma de matrizes facilitando o cálculo teórico.

Por último os gráficos de variação e conservação da energia total do sistema são gerados. 
Para o cálculo da conservação de energia total, é analisada a relação de temperatura com a variação de energia potencial.
As trajetórias registradas são usadas nos gráficos de movimento *Filmagem* e as **interações/transformações** são visualizadas para a variação de temperatura determinada pelo sistema.

