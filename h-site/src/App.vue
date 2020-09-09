<template>
	<div id="app">
		<div id="header">
			<Menu />
			<router-link to="/">Home</router-link>
			<SignLogIn />
		</div>
		<div id="nav">
			<router-view/>
		</div>
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
	margin-left: auto;
	margin-right: auto;
}
#header {
	background-color: black;
	width: 100%;
	margin: auto;
}

body {
	background-color: lightgrey
}

#nav {
	margin-left: auto;
	margin-right: auto;
	width: 500px;
}

#nav a {
	font-weight: bold;
	color: #2c3e50;
}

#nav a.router-link-exact-active {
	color: #42b983;
}
</style>
