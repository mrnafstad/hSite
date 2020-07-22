<template>
  <div id="app">
    <Menu />
    <div>
      <HelloWorld />
    </div>
  </div>
</template>

<script>
import HelloWorld from './components/HelloWorld.vue'
import Menu from './components/Menu.vue'
//import Logg from './components/Logg.vue'

export default {
  name: 'App',
  components: {
    HelloWorld,
    Menu
  }
}
</script>

<style>
#app {
  font-family: Avenir, Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  text-align: center;
  color: #2c3e50;
  margin-top: 60px;
  background-color: BlanchedAlmond;
  width: 400px;
  padding: 30px;
  border-radius: 10px;
  margin-left: auto;
  margin-right: auto;
}

  body {
    background-color: Beige;
  }

</style>
