<template>
  <div id="app">
    <Menu />
    <SignLogIn />
    <div id="nav">
      <router-link to="/">Home</router-link> 
    </div>
    <router-view/>
  </div>
</template>

<script>
import Menu from './components/Menu.vue'
import SignLogIn from './components/SignLogIn.vue'

export default {
  name: 'App',
  components: {
    Menu,
    SignLogIn
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
  background-color: lightgrey;
  width: 500px;
  padding: 30px;
  border-radius: 10px;
  margin-left: auto;
  margin-right: auto;
}

body {
  background-color: silver
}

#nav {
  padding: 30px;
}

#nav a {
  font-weight: bold;
  color: #2c3e50;
}

#nav a.router-link-exact-active {
  color: #42b983;
}
</style>
