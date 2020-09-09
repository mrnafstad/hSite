<template>
	<div id="header">
		<Menu />
		<router-link to="/">Home</router-link> | <router-link to="/blog">Blog</router-link>
		<SignLogIn />
	</div>
</template>

<script>
import Menu from './Menu.vue'
import SignLogIn from './SignLogIn.vue'

export default {
	name: 'MainHeader',
	components: {
		Menu,
		SignLogIn
	}
}
</script>

<style>
#header {
	background-color: black;
	width: 100%;
	margin: auto;
}

router-link {
	color: blue
}
</style>