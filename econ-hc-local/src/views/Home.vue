<template>
	<div>
		<div v-if="login">
			<h1>Velkommen {{ bruker.name }}</h1>
			<p>Utgifter: {{ bruker.utgifter }}</p>

			<p>Vil du lage en kobling med en annen bruker?</p>
			<ul>
				<li v-for="user in andreBrukere" :key="user.name"><button>{{ user.name }}</button></li>
			</ul>
			
		</div>

		<div v-else>
			<p>Ingen er logget inn</p>
			<button @click="$router.push('/login')">
				Logg inn på nytt</button>
		</div>
	</div>
</template>

<script>
	import store from '../store'

	export default {
		data() {
			return {
				bruker: {},
				andreBrukere: [],
				login: false
			}
		},
		async mounted() {
			//hentes etter at view endres. Må vente
			if (store.getters.loggin === true) {
				this.bruker = store.getters.bruker
				console.log(this.bruker.name)
			} else {
				this.bruker = await store.getters.bruker
				this.login = true
				this.andreBrukere = await store.getters.andre
				console.log(this.bruker.name)
				console.log(this.andreBrukere)
			}
		}
	}
</script>

<style>
	h1 {
		font-size: 15px;
	}
	li {
		/*background-color: silver;*/
		list-style-type: none;
		padding: 5px;
		display: inline-block;
		grid-template-columns: repeat(3, 1fr);
		grid-gap: 1rem;
		text-align: center;
		positon: relative;
		cursor: pointer;
	}

</style>
