<template>
	<div id="login-box">
		<h1>Login</h1>

		<p>Brukernavn:</p>
		<input id="input-box" v-model="brukernavn" placeholder="Brukernavn">

		<p>Passord:</p>
		<input class="input-box" type="password" v-model="passord" placeholder="*******">

		<div id="login-but">
			<button v-on:click="login">
				Logg inn
			</button>
		</div>
	</div>
</template>

<script>
//import firebase from 'firebase'
const db = require('../firebaseConfig.js')
import store from '../store'

	// // Your web app's Firebase configuration
	// var firebaseConfig = {
	// 	apiKey: "AIzaSyCg749u5HKVAVJ_tJqmnD1PWNNfTeQUBd4",
	// 	authDomain: "econ-hc.firebaseapp.com",
	// 	databaseURL: "https://econ-hc.firebaseio.com",
	// 	projectId: "econ-hc",
	// 	storageBucket: "econ-hc.appspot.com",
	// 	messagingSenderId: "762717017371",
	// 	appId: "1:762717017371:web:ebc9b937dd58e28efc2fbe"
	// };
	// // Initialize Firebase
	// firebase.initializeApp(firebaseConfig);

	// var db = firebase.firestore();
	//var users = db.collection("users");


export default {
	data() {
		return {
			navn: '',
			utgifter: 0,
			brukernavn: '',
			passord: ''
		}
	},
	methods: {
		login() {
			if (!store.getters.login) {
				db.usersCollection.where("brukernavn", "==", this.brukernavn).where("password", "==", this.passord).get().then(function(querySnapshot) {
					querySnapshot.forEach(function(doc) {
						store.commit('setPerson', doc.data())
						console.log(doc.id, " => ", doc.data())
					})
				}).catch(function(error) {
					console.log("Error getting documents: ", error)
				})
				//console.log(bruker.navn)
				this.hentAndreBrukere()
			}
			//sett brukernavn og passord, spør firebase om godkjenning
			this.$router.push('/kalk')
		},
		hentAndreBrukere() {
			db.usersCollection.get().then(function(querySnapshot) {
				querySnapshot.forEach(function(doc) {
					store.commit('setEksterneBrukere', doc.data())
					console.log(doc.id, " => ", doc.data())
				})
			}).catch(function(error) {
				console.log("Error getting documents: ", error)
			})
			console.log(store.getters.andre)
		}
	}
}
	
</script>

<style>
	#login-box {
		background-color: LightGrey;
		margin: 0;
		padding: 2em;
		position: absolute;
		border-radius: 1em;
		top: 50%;
		left: 50%;
		margin-right: -50%;
		transform: translate(-50%, -50%)
	}

	#login-but {
		padding-top: 1em;
	}
	button {
		border-radius: 0.25em;
		padding: 0.5em;
		background-color: red;
	}

	#link {
		border-radius: 0.15em;
		padding: 0.5em;
		color: black;
		background-color: red;
		text-decoration: none;
	}

	#input-box {
		border-radius: .25em;
	}

	h1 {
		font-size: 15px;
	}
</style>
